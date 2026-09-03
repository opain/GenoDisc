#!/usr/bin/env python3
"""Summarise Snakemake benchmark files by rule and flag rules relying on
profiles/slurm/config.yaml's default-resources (runtime/mem_mb) whose
observed usage leaves little or no cushion.

Rules that declare their own `resources:` block (get_mem_mb_fine,
mem_mb=20000, etc.) are skipped - they're sized individually already, not
via the shared default.

Usage:
    python scripts/analyze_benchmarks.py
    python scripts/analyze_benchmarks.py --jobs-dir /path/to/jobs --cushion 0.3
"""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

PIPELINE_DIR = Path(__file__).resolve().parent.parent
RULE_RE = re.compile(r'^rule (\w+):', re.MULTILINE)
BENCHMARK_RE = re.compile(r'benchmark:\s*\n\s*(f?"[^"]*"|f?\'[^\']*\')')
WILDCARD_RE = re.compile(r'\{(\w+)\}')

# Snakemake benchmark columns don't include disk usage - only wall time and
# memory are checked against default-resources here.


def parse_rules(pipeline_dir):
    """Returns (default_resource_rules, overridden_rules): {rule_name: benchmark_basename_template}."""
    default_rules = {}
    overridden_rules = {}

    for smk_file in sorted((pipeline_dir / 'rules').glob('*.smk')):
        text = smk_file.read_text()
        starts = [(m.group(1), m.start()) for m in RULE_RE.finditer(text)]
        for i, (rule_name, start) in enumerate(starts):
            end = starts[i + 1][1] if i + 1 < len(starts) else len(text)
            block = text[start:end]

            bm_match = BENCHMARK_RE.search(block)
            if not bm_match:
                continue  # no benchmark file -> nothing to analyze
            raw = bm_match.group(1)
            is_fstring = raw.startswith('f')
            literal = raw[2:-1] if is_fstring else raw[1:-1]
            if is_fstring:
                # f-string: {var} already interpolated by Python: only {{wildcard}}
                # survives as a literal Snakemake wildcard -> unescape it.
                literal = literal.replace('{{', '\x00').replace('}}', '\x01')
                literal = re.sub(r'\{\w+\}', '', literal)  # drop real f-string exprs (e.g. {resdir})
                literal = literal.replace('\x00', '{').replace('\x01', '}')
            basename = literal.rsplit('/', 1)[-1]

            has_resources = re.search(r'^\s*resources:', block, re.MULTILINE) is not None
            (overridden_rules if has_resources else default_rules)[rule_name] = basename

    return default_rules, overridden_rules


def compile_matcher(basename_template):
    pattern = re.escape(basename_template)
    pattern = re.sub(r'\\\{(\w+)\\\}', r'(.+)', pattern)
    return re.compile(f'^{pattern}$')


def collect_benchmarks(jobs_dir, rule_matchers):
    """Returns {rule_name: [(seconds, max_rss_mb), ...]}, unmatched: [Path, ...]."""
    usage = defaultdict(list)
    unmatched = []

    for tsv_path in sorted(jobs_dir.glob('*/benchmarks/*.tsv')):
        name = tsv_path.name
        matched_rule = None
        for rule_name, matcher in rule_matchers.items():
            if matcher.match(name):
                matched_rule = rule_name
                break
        if matched_rule is None:
            unmatched.append(tsv_path)
            continue

        try:
            with tsv_path.open() as f:
                row = next(csv.DictReader(f, delimiter='\t'))
            usage[matched_rule].append((float(row['s']), float(row['max_rss'])))
        except (StopIteration, ValueError, KeyError):
            continue  # empty/malformed benchmark file (e.g. job killed mid-write)

    return usage, unmatched


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pipeline-dir', type=Path, default=PIPELINE_DIR)
    ap.add_argument('--jobs-dir', type=Path, default=Path('/scratch/prj/neurohackpain/GenoDisc/jobs'))
    ap.add_argument('--runtime-default-min', type=float, default=480,
                     help='default-resources runtime in minutes (profiles/slurm/config.yaml)')
    ap.add_argument('--mem-default-mb', type=float, default=4000,
                     help='default-resources mem_mb (profiles/slurm/config.yaml)')
    ap.add_argument('--cushion', type=float, default=0.2,
                     help='fraction of the default a rule must stay under to count as having headroom (default 0.2 = 20%%)')
    ap.add_argument('--show-unmatched', action='store_true',
                     help='list benchmark files that could not be matched to any rule')
    args = ap.parse_args()

    default_rules, overridden_rules = parse_rules(args.pipeline_dir)
    all_matchers = {r: compile_matcher(t) for r, t in {**default_rules, **overridden_rules}.items()}
    usage, unmatched = collect_benchmarks(args.jobs_dir, all_matchers)

    runtime_limit_s = args.runtime_default_min * 60
    mem_limit_mb = args.mem_default_mb

    rows = []
    for rule_name in default_rules:
        samples = usage.get(rule_name)
        if not samples:
            continue
        times = [s for s, _ in samples]
        rss = [r for _, r in samples]
        max_s, mean_s = max(times), sum(times) / len(times)
        max_rss, mean_rss = max(rss), sum(rss) / len(rss)

        flags = []
        if max_s > runtime_limit_s:
            flags.append('RUNTIME EXCEEDED')
        elif max_s > runtime_limit_s * (1 - args.cushion):
            flags.append('runtime low cushion')
        if max_rss > mem_limit_mb:
            flags.append('MEM EXCEEDED')
        elif max_rss > mem_limit_mb * (1 - args.cushion):
            flags.append('mem low cushion')

        rows.append((rule_name, len(samples), max_s, mean_s, max_rss, mean_rss, flags))

    rows.sort(key=lambda r: r[2], reverse=True)  # by max runtime desc

    header = f"{'rule':<40}{'n':>4}{'max_min':>10}{'mean_min':>10}{'max_MB':>10}{'mean_MB':>10}  flags"
    print(header)
    print('-' * len(header))
    for rule_name, n, max_s, mean_s, max_rss, mean_rss, flags in rows:
        flag_str = ', '.join(flags)
        print(f"{rule_name:<40}{n:>4}{max_s/60:>10.1f}{mean_s/60:>10.1f}{max_rss:>10.0f}{mean_rss:>10.0f}  {flag_str}")

    flagged = [r for r in rows if r[6]]
    print(f"\n{len(rows)} default-resource rules with benchmark data "
          f"(runtime default={args.runtime_default_min:.0f} min, mem default={mem_limit_mb:.0f} MB, "
          f"cushion={args.cushion:.0%}).")
    print(f"{len(flagged)} rule(s) flagged: {', '.join(r[0] for r in flagged) or 'none'}")

    no_data = sorted(set(default_rules) - set(usage))
    if no_data:
        print(f"\n{len(no_data)} default-resource rule(s) with no benchmark data found under {args.jobs_dir}: "
              f"{', '.join(no_data)}")

    if args.show_unmatched and unmatched:
        print(f"\n{len(unmatched)} benchmark file(s) not matched to any rule:")
        for p in unmatched:
            print(f"  {p}")


if __name__ == '__main__':
    main()
