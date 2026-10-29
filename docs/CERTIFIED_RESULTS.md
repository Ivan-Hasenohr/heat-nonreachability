# Certified reference results

Reference campaign: `20261029_001624_932`, MATLAB R2023b Update 1 and INTLAB 12.
Six examples and 40 L2-table entries were certified. Total runtime was about 72 minutes.
These are recorded campaign results, not a new computation.

| Case | Nx | Nt | Lower J | Upper J |
|---|---:|---:|---:|---:|
| theorem6 | 800 | 500000 | -0.0010724613792912647 | -0.00023842757537720691 |
| corollary18 | 1200 | 3200000 | -0.00048891744109454008 | -0.00026764078305542511 |
| corollary19 | 1200 | 400000 | -0.00017983116379290289 | -0.00014744315584171668 |
| proposition20 | 1200 | 3200000 | -0.001897604780306196 | -0.0016612923852818444 |
| proposition22 | 500 | 25000 | -0.0094914346215048007 | -0.0034070759571060774 |
| proposition24 | 1000 | 80000 | -0.015902364821094694 | -0.011644186874956969 |

The decimal endpoints above reproduce the saved floating endpoints for inspection.
Use the original MAT interval objects for verified calculations; any rounded
intervals quoted in the article should be widened outwards.

Full budgets, table amplitudes and result hashes are in
[../data/certification_summary.json](../data/certification_summary.json).
For mode 1 the best table value is 0.1435 at Nx=1000, whereas Nx=2000 gives 0.1440.
Both are certified at their respective meshes; the table is not monotone in this mode.

Large MAT files are distributed separately; see [RESULTS.md](RESULTS.md).
