# GenoDisc dataset contribution — template

Fill this in when contributing a new **reference dataset** to the GenoDisc pipeline (an LD reference,
annotation, weights panel, gene map, etc. staged into the pipeline's resource dir). Copy this file to
`pipeline/docs/datasets/<name>.md` and complete every field. Leave no `TODO`s.

The contributor adds a `download_<name>` rule (see `CONTRIBUTING.md` → "Adding a reference dataset");
the maintainer stages it on the deployment's resources and, if the dataset is user-selectable, wires
it into the website. See `rules/dependencies.smk` (`download_ldsc_scores`, `download_eqtlgen`,
`download_biomart`) for the staging pattern.

---

## 1. Dataset identity
- **Display name:** <!-- short human name -->
- **One-line summary:** <!-- what the data is and which method(s)/analysis consume it -->
- **Used by:** <!-- the method(s)/rule(s) that need it, e.g. sbayesrc -->

## 2. Source & provenance
- **Source / download URL:** <!-- stable OFFICIAL URL (prefer the primary host over mirrors/drives) -->
- **Version / release:** <!-- dataset version or release date -->
- **Genome build:** <!-- GRCh37 / GRCh38 -->
- **Ancestry / population(s):** <!-- EUR/EAS/AFR/... or N/A; one entry per published ancestry -->
- **Size:** <!-- download size + unpacked size (matters for staging/disk) -->
- **Checksum (if published):** <!-- md5/sha, else "none published" -->

## 3. Licence & redistribution
- **Licence:** <!-- data licence -->
- **Redistribution terms:** <!-- may it be mirrored/hosted? FLAG restrictions, e.g. UKB-derived data
  (SBayesRC UKB LD references) have usage restrictions — must download from source, not re-host -->
- **Citation users should give:** <!-- paper(s) to cite when using this data -->

## 4. Staging (pipeline)
- **Download rule:** `download_<name>` in `rules/dependencies.smk` — summarise what it does
  (wget/unzip/tar; or an R prep script). Prefer official URLs.
- **On-disk layout:** `{resdir}/data/<name>/...` <!-- the dir/files produced; note <POP>/ subdirs etc. -->
- **Sentinel output:** <!-- the concrete file or `.done` marker the rule declares as output -->
- **Config key:** <!-- T/F toggle, or an explicit path key (with NA = use the staged copy in resdir), or
  "none (core/always-staged)" -->
- **resource_inputs / prepare_resources:** <!-- which gated block it is appended to -->

## 5. User-selectable? (website)
- **Panel-picker needed:** <!-- yes if the user chooses it (like a TWAS/SMR panel) → maintainer adds a
  picker entry + docs row; no if it is an internal reference used by a method -->

## 6. Maintainer / contact & notes
- **Maintainer / contact:** <!-- name + email -->
- **Notes / caveats:** <!-- known issues, ancestry gaps, build-conversion needs, etc. -->
