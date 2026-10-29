# Computer-assisted non-reachability certificates

MATLAB/INTLAB code accompanying **Computer-assisted proofs of non-reachability for linear parabolic PDEs under bounded control constraints**, by Ivan Hasenohr, Camille Pouchol, Yannick Privat and Christophe Zhang.

Companion preprint: [HAL hal-05357694](https://hal.science/hal-05357694). The scripts use the revised-manuscript numbering associated with the certification campaign of 29 October 2026; earlier preprint versions may use different numbers.

The code separates three tasks:

1. **Floating-point search:** find a candidate dual witness using Chambolle–Pock.
2. **Verified evaluation:** evaluate a fixed witness using INTLAB and rigorous discretization bounds.
3. **Visualization:** illustrate the states, targets, control regions and witnesses.

A successful floating-point search is not a certificate. Non-reachability is certified only by the verified evaluation described below.

## Requirements and setup

The development checks used:

- MATLAB R2023b Update 1;
- INTLAB version 12, installed separately.

Other versions have not been systematically tested. INTLAB is not bundled with this repository.

Initialize INTLAB according to its installation instructions. Then change to the repository root and add the source directory to the MATLAB path. For example, replacing both paths:

```matlab
addpath('C:/path/to/INTLAB');
startintlab;

% Use MATLAB's native cleanup object, not INTLAB's legacy compatibility stub.
cleanup_file = which('onCleanup');
if contains(cleanup_file, [filesep 'newfun' filesep 'onCleanup' filesep])
    rmpath(fileparts(cleanup_file));
end

cd('C:/path/to/paper-cap');
project_root = pwd;
addpath(genpath(fullfile(project_root,'src')));
```

Add only the INTLAB root and call `startintlab`; avoid recursively adding the entire INTLAB installation. If its legacy `newfun/onCleanup` stub has already caused a test-framework error, use a fresh MATLAB session after correcting the path.

Avoid adding an older version of the project to the same MATLAB path. The functions used can be checked with:

```matlab
which init_heat_fem
which certify_heat_nonreachability
which certify_coupled_heat_nonreachability
```

## A small example

The following example exercises the scalar certificate evaluator on a small mesh:

```matlab
Nx = 20;
Nt = 100;

args = init_heat_fem(Nx,Nt,intval('1'),false);

x = (1:Nx-1)'/Nx;
p = sin(pi*x);
p = p/(norm(p)*sqrt(1/Nx));

y_0 = zeros(Nx-1,1);
y_f_base = intval('0.035');

[certified,J,~,~,~,bounds] = ...
    certify_heat_nonreachability(p,y_0,y_f_base,args);

disp(certified);
disp(J);
disp(bounds);
```

Here `y_f_base` defines the continuous target
\(y_f(x)=0.035\sin(\pi x)\).

This small mesh is intended for a quick execution check. It may be too coarse to certify non-reachability. Use the paper scripts below to reproduce the intended configurations.

## Reproducing the paper examples

Paths below are relative to the repository root. The numbering follows the revised manuscript.

| Result | Script | Result file |
|---|---|---|
| Theorem 6 | `examples/heat_positive/run_heat_intro_certificate.m` | `results/theorem6.mat` |
| Corollary 18 | `examples/heat_positive/run_heat_positive_certificate.m` | `results/corollary18.mat` |
| Corollary 19 | `examples/heat_min_time/run_heat_min_time_certificate.m` | `results/corollary19.mat` |
| Proposition 20 | `examples/heat_tent/run_heat_tent_certificate.m` | `results/proposition20.mat` |
| Proposition 22 | `examples/coupled_ball/run_coupled_ball_certificate.m` | `results/proposition22.mat` |
| Proposition 24 | `examples/coupled_partial/run_coupled_partial.m` | `results/proposition24.mat` |
| L2 table | `examples/table_l2/compute_heat_l2_table.m` | `results/table_l2/table_L2_all.mat` |

For example:

```matlab
run(fullfile(project_root,'examples','coupled_ball', ...
    'run_coupled_ball_certificate.m'));
```

The scripts contain the final time, mesh sizes, continuous data and control parameters. Run examples individually in a clean workspace.

Certification scripts use saved witnesses from `data/witnesses`, except where the witness is constructed directly. Reproducing a certificate does not require rerunning its optimization.

Scripts can overwrite their result files and transferred witness files. Preserve previous experiments before rerunning them. Fine-mesh computations can require substantial time and memory.

## What constitutes a certificate?

For each fixed witness, the verifier computes an interval enclosure of the continuous separation functional.

The returned flag is true only when:

```matlab
isfinite(inf(J)) && isfinite(sup(J)) && sup(J) < 0
```

A strictly negative upper endpoint proves non-reachability under the assumptions and estimates established in the article.

If `certified` is false, the evaluation is inconclusive. It does not establish reachability. Possible causes include an unsuitable witness, insufficient resolution or overly wide interval enclosures.

The witness need not be an optimizer. Floating-point search and interpolation produce candidate coefficients; the verified evaluation treats the resulting coefficients as fixed witness data.

## Saved results and error information

The six certificate scripts save the witness, problem arguments, continuous state/target data, certification flag, final interval and computation time.

The `bounds` structure contains:

| Field | Meaning |
|---|---|
| `time` | Upper bound on the temporal contribution to the discretization error |
| `space` | Upper bound on the spatial contribution |
| `discretization` | Combined discretization-error bound actually used |
| `discrete_enclosure` | Interval evaluation before adding discretization uncertainty |
| `roundoff` | Radius of `discrete_enclosure` |
| `final_enclosure` | Enclosure after adding discretization uncertainty |

`roundoff` is a diagnostic name: this radius can include interval overestimation and input uncertainty as well as floating-point rounding effects.

The combined bound is evaluated directly in interval arithmetic. Do not reconstruct it by adding separately printed, rounded contributions.

In the L2 table:

- `lambda_num_all` contains floating representations of the reported amplitudes;
- `lambda_cert_all` contains the interval representations used for their certification;
- `J_all` contains the final functional enclosures;
- `sups_all` contains the upper bounds used for the certification checks;
- `completed_rows` identifies the rows completed and saved.

The table is saved after each complete row. Automatic restart from an interrupted computation is not implemented.

## Data representations

Continuous data and floating nodal samples have different roles.

- `y_0` and `y_f` in search workflows contain floating nodal samples.
- `y_0_base` and `y_f_base` define the continuous data used by the verifier.
- `p` contains interior nodal witness coefficients. The verifier constructs its continuous cubic-spline witness from these values.

The current representation convention is selected by vector length:

| Representation | Scalar data length | Two-component data length |
|---|---:|---:|
| P1 nodal values | `Nx-1` | `2*(Nx-1)` |
| Cubic Hermite coefficients | `2*Nx` | `4*Nx` |
| Sine-series coefficients | Other supported length | Other supported even length |

For cubic Hermite data, each component contains interior function values followed by **mesh size times the derivatives**, including endpoint derivatives. Function values at the two endpoints are zero.

For sine-series data, coefficients represent
\[
f(x)=\sum_{k=1}^{m}a_k\sin(k\pi x).
\]

Coupled data concatenate the first component and then the second. Both components use the same representation and coefficient count.

Because selection is based on length, sine-series lengths can collide with P1 or Hermite lengths after changing `Nx`. The example scripts include assertions against these collisions.

For the tent target, `Nx` must also be divisible by 5 so its vertex is a mesh node.

## Witness normalization

The operation

```matlab
p = p/(norm(p)*sqrt(1/Nx));
```

uses the nodal norm. It does not impose an exact unit L2 norm on the continuous spline.

The original separation functional is positively homogeneous, so positive rescaling preserves its sign. Norms needed for the verified estimates and the target-ball contribution are evaluated separately.

For coupled witnesses, normalization is joint across both components. Their relative amplitudes must be preserved.

## Searching for new witnesses

Floating search scripts are provided alongside the certificate scripts:

- `examples/heat_positive/run_heat_positive.m`
- `examples/heat_min_time/run_heat_min_time_search.m`
- `examples/heat_tent/run_heat_tent_search.m`
- `examples/coupled_ball/run_coupled_ball_search.m`

Floating propagation tables are disabled by default:

```matlab
args = init_heat_fem(Nx,Nt,intval_T,false);
```

Searches require the floating table:

```matlab
[args,args_float] = init_heat_fem(Nx,Nt,intval_T,true);
```

The coupled initializer follows the same convention. Use coarse meshes for search; full floating propagation tables can be expensive.

The scalar workflow supports only a zero initial state. The coupled workflow supports nonzero initial states, partial targets and L2-ball targets. A ball radius must be passed explicitly to the coupled workflow and verifier.

Search diagnostics are not certificates. In particular, the coupled stopping diagnostic integrates the support of a reconstructed, masked P1 field while retaining nodal target terms. It can differ from the nodal functional associated with the search operators.

Regularization, when enabled in the scalar search, helps construct a witness. Certification evaluates the original, unregularized functional.

For the tent example, `examples/heat_tent/compare_heat_tent_witnesses.m` transfers the raw witness to the regularized result mesh, evaluates the original functional and saves both second-derivative norms and enclosures in `results/tent_witness_comparison.mat`. This is an additional fine-mesh evaluation; it is not part of the seven-job reference campaign.

## Source organization

- `src/setup`: initializers and argument containers;
- `src/search`: floating optimization routines and helpers;
- `src/workflows`: search followed by verified evaluation;
- `src/numerics`: integration, spline and propagation kernels;
- `src/verification`: verified supports, error bounds and certificate evaluators;
- `examples`: configurations used in the article;
- `data/witnesses`: saved candidate witnesses;
- `results`: saved evaluations;
- `tests`: regression tests and reference routines;
- `figures`: Python figure generator and portable floating snapshots.

The contents and types of `args` are documented in `src/setup/make_control_args.m`. This function only stores arguments; it does not validate or certify arbitrary inputs.

Verified arguments must come from the appropriate initializer. Do not substitute `args_float` for `args`, or modify related cells independently without checking their consistency.

Propagation routines evaluate implicit-Euler powers, not exact matrix exponentials. Their discrepancy from continuous evolution is covered by separate discretization estimates.

## CAP regression tests

After setup, run:

```matlab
addpath(fullfile(project_root,'tests'));
test_results = run_cap_tests;
```

The suite contains 30 grouped tests covering numerical kernels, control supports,
error bounds, certificate decisions and input contracts. See
[tests/README_CAP.md](tests/README_CAP.md) for the independent references and limits.
These small tests do not overwrite paper results or replace the mathematical estimates.

## Complete certification campaign

```matlab
addpath(fullfile(project_root,'examples'));
campaign_dir = run_cap_campaign(true);
```

This runs the six certificates and then the complete L2 table. Pass `false` to omit
the table. The launcher copies sources and input witnesses before computing.
New results are saved under `results/campaigns/<timestamp>/code/results`, with
`summary.mat` and per-case `logs/` beside `code/`. An error in one case is recorded
and does not prevent the remaining cases from running. Completed campaigns are not
resumed automatically.

Standalone scripts instead write directly under `results/`. A campaign does not
replace those root result files. The helper `examples/check_campaign_result.m`
checks saved formats and flags; it does not recompute a proof.

## Reference results and figures

[docs/CERTIFIED_RESULTS.md](docs/CERTIFIED_RESULTS.md) records the successful
29 October 2026 campaign. The small summary and floating plotting snapshots are
included in Git; large MAT files and campaign copies are excluded.
See [docs/RESULTS.md](docs/RESULTS.md) for the separate results archive and its layout.

The paper figures are generated with Python:

```console
python -m pip install -r figures/requirements.txt
python figures/run_paper_figures.py
python figures/test_paper_figures.py
```

Default figure generation uses bundled floating snapshots of the reference campaign
and needs neither MATLAB nor INTLAB. Exports go to `results/figures/`.
See [figures/README.md](figures/README.md) for refreshing data, display scales and
the optional notebook. Drawn curves and floating snapshots are illustrations,
not verified error enclosures.

## License and citation

The source code is licensed under the **PolyForm Noncommercial License 1.0.0**;
see [LICENSE](LICENSE). Use, modification and redistribution are allowed for the
purposes permitted by that license. Copies must retain the supplied `Required
Notice`, which identifies the companion article and its authors.

The license explicitly permits uses by educational institutions and public research
organizations regardless of their funding source. Commercial uses outside the
license's permissions require separate authorization; consult the actual license
text rather than treating "noncommercial" as a restriction on every funded project.

For publications using the code, please cite the companion article as specified in
[CITATION.cff](CITATION.cff). This academic citation is requested; the binding notice
requirement concerns redistribution and is defined in LICENSE.

MATLAB and INTLAB are external dependencies and are not redistributed here.
