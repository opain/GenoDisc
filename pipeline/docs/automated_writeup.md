# Automatically generated draft write-up

Every GenoDisc run now produces a **draft methods section, a descriptive results
section and a matching reference list**, generated from what the pipeline actually
ran. They are written into the results bundle alongside the data the Output Viewer
reads.

> **These are drafts, not finished text.** They are a starting point to save you
> re-deriving tool versions, parameters, test counts and significant-result counts
> by hand. You must read, check and edit them before using them in a paper. Every
> generated file opens with a banner saying exactly this.

## Where the files are

They can be downloaded from the **References tab of the Output Viewer** ("Draft
write-up" panel: Methods, References, and an All-files `.tar.gz`), and they also live
inside the results bundle (`bundle.tar.gz`) under `package/writeup/`:

| File | Contents |
|---|---|
| `methods.md` | One paragraph per analysis that ran, in a fixed order, naming the tool, its version, the reference data and the key parameters, plus a QC/harmonisation paragraph and a closing pipeline/web-server paragraph. |
| `<GWAS>/results.md` | A strictly descriptive results section for each submitted GWAS: number of tests, number significant, the top result and a table/tab reference, with explicit sentences for null results. |
| `references.bib` | BibTeX entries for exactly the tools and datasets this run used. |
| `references.ris` | The same references in RIS format. |

The `manifest.json` records the generated filenames, the generator version and the
results-schema version it read, under a `writeup` key.

## What the numbers mean

- **"Significant"** always means the stored Benjamini–Hochberg FDR-adjusted
  P-value (`P.FDR`) is below 0.05 — the same column and threshold the Output Viewer
  applies by default, so the counts reconcile with the viewer. For TWAS/PWAS/SMR the
  text additionally reports how many of the FDR-significant features were
  *colocalised* (FUSION) or *HEIDI-supported* (SMR); these stricter counts match the
  viewer's "high-confidence" Overview tiles.
- **Fine-mapping** has no P-value: "significant" means a 95% credible set was
  resolved.
- **Heritability** (LDSC / SBayesRC / SumHer) has no significance test; the text
  reports the estimate and its uncertainty. The three methods use different models
  and are not expected to agree.
- Analyses that were **switched off** are omitted entirely (never described as "not
  performed"). Analyses that were **requested but produced no output** for any GWAS
  are listed in a *Reconciliation — check before use* footer, with the likely reason
  (e.g. a European-only step run on a non-European GWAS, or a pipeline failure to
  check in the logs).

## Known limitations (read before publishing)

- **Version gaps.** The MAGMA and SMR binaries carry no version string anywhere
  (no version in the download URL, no checksum), and `plink2R` is pinned to a moving
  `master`. Where a version cannot be recovered, the text says *"version not pinned
  upstream"* rather than inventing one — fill these in from your records if needed.
- **`[CITE: GenoDisc paper]` placeholders.** The closing paragraph cites the GenoDisc
  pipeline and web server with placeholders; replace them with the real citation(s).
- **Citations to verify.** A few entries in `references.bib` are marked `% VERIFY`
  (they were not in the maintainer-curated in-app reference table); confirm them.
- **No interpretation.** The results text is deliberately descriptive — it reports
  counts and top hits only, with no biological interpretation, causal or novelty
  claims. Add interpretation yourself.
- **Drafts can drift from intent.** The generator describes what *ran*, which may
  differ from what you meant to run; the reconciliation footer is there to catch
  that.

## For maintainers

The generator is pipeline-side (`scripts/functions/writeup/`), invoked from
`scripts/package_results.R` at packaging time. Its citation/version source of truth
is `writeup/references.bib` + `writeup/registry.yaml`, and fixed algorithm constants
live in `writeup/params.yaml`. When you add a new method or dataset, add a
`registry.yaml` entry (with a `bibkey`) and the BibTeX entry, and — if it introduces
a new block — a spec in `writeup.R` plus a methods/results builder. The tests
(`tests/testthat/test-writeup-*.R`) enforce that every block-enabling config flag
maps to a registry entry with a real `references.bib` key, and snapshot the full
fixture output so template changes are visible in diffs.
