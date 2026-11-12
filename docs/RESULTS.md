# Campaign outputs

The cleaned repository starts with fixed witness inputs and no numerical
certification outputs. Generate a new reference set with `run_cap_campaign(true)`.

Each dated directory under `results/campaigns/` contains:

- `initial_witnesses/`, the five original inputs before running any script;
- `code/`, the frozen sources, examples, witness inputs and their result files;
- `logs/`, one transcript for each of seven certificates and the L2 table;
- `summary.mat`, case statuses, runtime information and output paths.

`summary.all_certified` is true only if every selected job succeeded.
The MAT files contain fixed witnesses, target data and verified intervals.
The table file records completed rows and the certified decimal amplitudes.
An interrupted campaign keeps its completed results for inspection.

After success, keep the full dated campaign as the reference for the paper.
Generate figures from its `code/` directory as documented in
[PAPER_FIGURES.md](PAPER_FIGURES.md). Historical results are in the local
archive created before cleanup.
