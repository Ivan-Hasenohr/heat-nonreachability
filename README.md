# Computer-assisted non-reachability certificates

Code accompanying *Computer-assisted proofs of non-reachability for linear parabolic PDEs under bounded control constraints*, by Ivan Hasenohr, Camille Pouchol, Yannick Privat and Christophe Zhang.

Companion preprint: [HAL hal-05357694](https://hal.science/hal-05357694).
The examples below follow the resubmitted manuscript. Filenames describe the
physical cases so they remain usable if the paper numbering changes.

## Requirements

The verification environment is MATLAB R2023b Update 1 and INTLAB 12.
Initialize your separately installed INTLAB, then add the project sources:

```matlab
addpath('C:/path/to/INTLAB');
startintlab;
cd('C:/path/to/paper-cap');
project_root = pwd;
addpath(genpath(fullfile(project_root,'src')));
```

Add only the INTLAB root before `startintlab`. Avoid adding another copy of
paper-cap to the same MATLAB path. The test runner handles INTLAB's legacy
`onCleanup` compatibility path if it masks MATLAB's native implementation.
MATLAB and INTLAB are external dependencies.

## Paper examples

| Result | Physical case | Script in examples/ | Standalone result in results/ |
|---|---|---|---|
| Theorem 6 | heat_intro | heat_positive/run_heat_intro_certificate.m | heat_intro.mat |
| Corollary 20 | heat_positive | heat_positive/run_heat_positive_certificate.m | heat_positive.mat |
| Corollary 21 | heat_min_time | heat_min_time/run_heat_min_time_certificate.m | heat_min_time.mat |
| Proposition 22 | heat_projected | heat_projected_target/run_heat_projected_target_certificate.m | heat_projected.mat |
| Proposition 24 | coupled_ball | coupled_projected/run_coupled_projected_ball.m | coupled_projected/ball/fine.mat |
| Proposition 26 | coupled_partial | coupled_projected/run_coupled_projected_partial.m | coupled_projected/partial/fine.mat |
| Proposition 27 | coupled_split | coupled_projected/run_coupled_projected_split.m | coupled_projected/split/fine.mat |
| Table 1 | table_l2 | table_l2/compute_heat_l2_table.m | table_l2/table_L2_all.mat |

For example:

```matlab
run(fullfile(project_root,'examples','coupled_projected','run_coupled_projected_split.m'));
```

Each script fixes the model, exact target, witness and mesh. The coupled point
target uses radius zero. The partial target has a free second component and
a witness whose second component is exactly zero.

## Complete campaign

```matlab
addpath(fullfile(project_root,'examples'));
campaign_dir = run_cap_campaign(true);
```

This runs the seven certificates and all 40 table entries. `false` skips the
table. Before computing, the launcher freezes the sources, example scripts
and five witness inputs under `results/campaigns/<timestamp>/code/`.
It saves a `summary.mat` and one log per case. A failed case is recorded and
the remaining cases continue. There is no automatic resume of an interrupted
campaign. Fine meshes require substantial time and memory.

Generate all seven figures from that same completed campaign:

```console
python -m pip install -r figures/requirements.txt
python figures/run_paper_figures.py --project "results/campaigns/<timestamp>/code"
```

The generated PDF, PNG and SVG files and `figure_metadata.json` are placed
inside that campaign's `code/results/figures/`. Every plotted witness is
the natural cubic spline reconstructed from the saved nodal values.

## Fixed inputs and search

`data/witnesses/` contains five fixed inputs:

- `heat_intro_seed.mat`, the nodal seed transferred to Nx=800 for Theorem 6;
- `heat_positive_coarse.mat`, the coarse witness for Corollary 20;
- `heat_min_time_coarse.mat`, the coarse witness for Corollary 21;
- `heat_projected_sine.mat`, the 24-mode scalar seed for Proposition 22;
- `coupled_split_sine.csv`, the 16-mode two-component seed for Proposition 27.

The ball and partial-target seeds are stated directly in their shared launcher.
Reproducing the certificates does not require rerunning witness searches.
The two optional scalar search scripts save their experiments under
`results/search/`; fixed input witnesses are kept separate from those outputs.

## Verification and tests

A certificate succeeds when its final interval has finite endpoints and a
strictly negative upper endpoint. An unsuccessful sign check is inconclusive.
The witness need not minimize the separation functional.

```matlab
addpath(fullfile(project_root,'tests'));
test_results = run_cap_tests;
test_heat_sine_filters;
```

The MATLAB tests use small meshes and independent mathematical references.
Scalar comparison and norm checks can also be run with
`examples/heat_checks/run_heat_sine_checks.m`.

```console
python figures/test_paper_figures.py
```

The plotting tests use synthetic inputs and independent integration formulas,
so they can run before the first complete campaign.

## Numerical data and bounds

`bounds.time` and `bounds.space` record the discretization contributions.
`bounds.discretization` is the combined bound actually used in the proof.
`bounds.discrete_enclosure` and `bounds.final_enclosure` enclose the discrete
and continuous functionals. `bounds.roundoff` is the radius of the discrete
enclosure, including input uncertainty and interval overestimation.
Use saved INTLAB objects for computations and widen decimal intervals
outwards when reporting them. Nodal normalization does not impose an exact
unit continuous L2 norm on the spline.

Continuous state/target data use the existing length-based representations:
P1 nodal values, cubic Hermite coefficients or sine coefficients. Example
assertions reject collisions of their lengths. The propagation routines
evaluate implicit-Euler powers, with separate bounds for the continuous PDE.

Before cleanup, the complete previous repository, including results and Git
metadata, was archived locally. Fresh campaign outputs are kept under
`results/` and are excluded from Git. See [docs/RESULTS.md](docs/RESULTS.md).

## License and citation

The source uses the PolyForm Noncommercial License 1.0.0. Retain the required
notice in [LICENSE](LICENSE) when redistributing it. For publications using
the code, cite the article described in [CITATION.cff](CITATION.cff).
