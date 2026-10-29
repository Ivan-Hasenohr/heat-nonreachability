# Small CAP regression suite

These tests exercise numerical kernels, error budgets, certificate decisions and input contracts. They use small meshes and fixed witnesses. They do not optimize a witness, run the article's full examples, or overwrite saved witnesses or paper results.

## Run

MATLAB R2023b and INTLAB 12 are the tested versions. MATLAB's built-in `matlab.unittest` framework is used; no additional testing toolbox, Python, external fixture file, or audit workspace is required.

Start your own INTLAB installation, then run:

```matlab
addpath('C:/path/to/INTLAB');
startintlab;
addpath('C:/path/to/paper-cap/tests');
results = run_cap_tests;
```

The current working directory does not matter. The runner finds `src` relative to its own file and restores the MATLAB path afterwards. On failure, MATLAB reports the test name, failed assertion and source line; the runner raises an error. A single test file can also be run with `runtests` after adding `src` recursively to the path.

The suite creates one temporary MAT file to check save/load and removes it afterwards. It does not write a permanent test report. To keep results, save `results` to a location of your choice.

The intentionally unsupported one-cell control geometry produces warnings before its expected rejection. These warnings are normal in that test.

## Read and maintain

Each file groups one part of the proof computation. Each local `test...` function is an independently reported test. Loops cover multiple values of the same mathematical scenario; the number of reported tests is therefore smaller than the number of configurations checked.

| File in `tests/cap` | Purpose and independent references |
|---|---|
| `test_cap_numerics.m` | Euler powers, state/adjoint distinction, polynomial integrals, natural spline and sine pairings. References are closed-form powers and rational integrals. |
| `test_cap_control_supports.m` | Aligned/partial boundary cells, sign changes, multiple control intervals, interval witnesses, controlled component and block transitions at 2048/8192 columns. A separate cell-by-cell oracle integrates affine fields; geometric time factors and a nilpotent coupling test propagation. |
| `test_cap_error_bounds.m` | Constants C3, control measure, actual box amplitude, target projection error, transposed reaction, separate regularity denominators and optional diagnostics. Polynomial witnesses allow explicit norms of A* p. |
| `test_cap_certificates.m` | P1/spline/sine targets, both coupled target types, continuous product norm for ball radii, nonzero initial states, uncertainty inflation, decisions and metadata save/load. Initial propagation is recomputed from the weak FE system. |
| `test_cap_input_contracts.m` | Floating/verified separation, step consistency, matrix size, malformed witnesses/states, exact scalar zero initial state, explicit ball radius and optional floating tables. |

The four small helpers in `private/` only compare references, check expected rejections and construct explicit polynomial/spline fixtures. They do not call a production certifier to decide the expected sign.

Start with the numerical tests, then supports, error bounds, certificates and input contracts. A failure is a reason to inspect a specific calculation or contract, not evidence by itself that a theorem is false. If intentionally changing a formula or interface, revise its mathematical reference and test together. Do not relax a tolerance merely to obtain a pass.

## Reference calculations

For `q(x)=x(1-x)`, the suite uses

```text
integral(q^2) = 1/30,   ||q''|| = 2.
```

With the P1 hat centered at 1/2 and half-width 1/4, its squared norm is 1/6 and its pairing with q is 23/384. Its sine pairing follows by integrating the two affine branches.

For the natural spline through the nodal hat on Nx=8, the h-scaled derivatives are

```text
[-3, 6, -21, 78, 0, -78, 21, -6, 3] / 97.
```

They can be obtained independently by imposing zero endpoint second derivatives and solving

```text
m_(j-1) + 4*m_j + m_(j+1)
    = 6*(v_(j+1) - 2*v_j + v_(j-1))/h^2.
```

Piecewise polynomial integration then gives spline norm squared `71989/658630`, second-derivative norm squared `712704/97`, and P1-hat/spline pairing `731/7760`. The fixture stores these rational values rather than reading outputs previously saved by the production code.

For a coupled polynomial `p=u*q`, let `r=R'*u` and `k=diag(kappa)*u`. Direct integration gives

```text
||A* p||^2 = 4*sum(k.^2) - (2/3)*sum(k.*r) + sum(r.^2)/30.
```

This checks the direction of the reaction transpose independently of the implementation's Gram matrices. Some error-bound fixtures change reaction/form constants without rebuilding A: they test the bound helper alone and are not presented as complete PDE configurations.

## What passing means

- Positive final values, zero values and non-finite upper bounds must never be reported as certificates.
- A negative discrete functional can remain inconclusive after error inflation; this situation is tested deliberately.
- `bounds.roundoff` is the radius before discretization inflation, including possible interval overestimation and input uncertainty. It is not a measurement of rounding error alone.
- Narrow interval references are checked for overlap and numerical agreement. This is not a formal proof of full containment of one interval by another.
- For mixed-sign coupled support cells, only the upper bound is required to be tight; a deliberately conservative lower bound is allowed.
- Samples from interval witnesses check several particular fields, not every field in the interval box.
- Synthetic A=0, A=-I and nilpotent fixtures isolate kernel properties. They do not claim to be heat-equation proof configurations.
- The self-adjoint branch test also checks its explicit diagnostic label; changing that label requires updating this test.

This readable subset complements the broader audit suite. It does not include the full L2-table script/checkpoint harness, exhaustive geometry tests or workflow/memory-guard checks from that audit. It does not prove the correctness of arbitrary MATLAB code or replace the mathematical estimates, INTLAB's guarantees and the final six certificate runs plus full L2 table.

Before release, run this suite on the frozen sources, then reproduce the paper's results using their documented certificate scripts. Changing control parameters intentionally may require adapting tests that exercise the example initializers' current defaults.
