# Reference results

The reference certification campaign is `20261029_001624_932` (29 October 2026).
It certified six examples and all 40 entries of the L2 table.
See [CERTIFIED_RESULTS.md](CERTIFIED_RESULTS.md) and
[../data/certification_summary.json](../data/certification_summary.json).

## Separate archive

Large results are deliberately excluded from Git. The prepared distribution asset
is named `certified-results-20261029.zip`. Obtain that asset from the maintainers
until it is attached to a repository release; no release download URL has yet been
configured in this repository.

Extract the archive at the repository root. It provides:

- `results/theorem6.mat`, both corollaries and the three proposition MAT files;
- `results/table_l2/table_L2_all.mat`;
- `results/provenance/summary.mat`, original logs, source manifest and the frozen
  campaign sources under `results/provenance/code/`.

The MAT files contain the fixed witnesses, verified arguments and functional
enclosures; the six certificates also contain `bounds`. Loading INTLAB objects
requires INTLAB. The provenance copy is historical reference material; initialize
the active project with `src/`, not by recursively adding the whole repository.

The maintainers' copy of the distribution asset is under `results/distribution/`.
This directory and extracted results remain excluded from Git.

## Regenerate instead

After the setup described in the root README:

```matlab
addpath(fullfile(project_root,'examples'));
campaign_dir = run_cap_campaign(true);
```

The new campaign stores its own results under
`results/campaigns/<timestamp>/code/results`; it does not update root result files.
The reference campaign took about 72 minutes on the original machine, not a runtime
guarantee for other systems.
