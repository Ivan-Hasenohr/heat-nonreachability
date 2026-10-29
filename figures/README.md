# Paper figures in Python

This is the approved Matplotlib presentation of the five paper examples, refreshed
from campaign `20261029_001624_932`. Floating snapshots are bundled so that plotting
requires neither MATLAB nor INTLAB. They do not contain a CAP proof.

## Run from the repository root

Python 3.10 or newer:

```console
python -m pip install -r figures/requirements.txt
python figures/run_paper_figures.py
python figures/test_paper_figures.py
```

Five figures are exported as PDF, SVG and PNG under `results/figures/`. Existing
exports with the same names are replaced. `figure_metadata.json` records source
hashes and display factors; `captions.tex` supplies matching LaTeX captions.
The notebook is optional; use a Jupyter installation to open `paper_figures.ipynb`.
It contains no embedded copies of older figures and regenerates the current ones.

## Refresh from saved MATLAB results

After extracting the separate results archive at the repository root:

```console
python figures/run_paper_figures.py --project .
```

Or use the `code` directory of a completed campaign:

```console
python figures/run_paper_figures.py --project "results/campaigns/<timestamp>/code"
```

Refresh reads five result MAT files plus
`data/witnesses/proposition20_unregularized.mat`. It updates the floating NPZ
snapshots in `figures/data/`, then regenerates the figures. It does not change any
MATLAB source, witness MAT or certified result MAT. `--output` and `--data-dir`
can select different export and snapshot directories.

## Interpretation and customization

- Corollaries 18 and 19 have distinct witnesses and times. Initial state is green,
  target red, witness blue, maximal-control response black and control regions orange.
- Proposition 20 compares the raw search witness and the regularized witness at
  the same displayed continuous L2 norm. The raw search file contains no certificate;
  this does not mean that its transferred witness cannot be certified. A subsequent
  fine-mesh evaluation has also produced a negative enclosure.
- Coupled vectors use one common factor for their two components. Control acts on
  the second component. The proposition-22 target is an L2 ball, not a pointwise band.
- Proposition 24 has a prescribed first target and a free second target component;
  the latter is not represented as a zero target.

Display scaling never changes saved witnesses. Exact factors are recorded in
metadata and captions, including the separate enlargement of the uncontrolled
final state in proposition 22. The continuous spline norm used for display is
computed in floating arithmetic, not as a verified interval norm.

The maximal-control response uses a piecewise-quadratic stationary solution minus
a Fourier evolution; the coupled free evolution uses 2-by-2 modal exponentials.
These reference curves and all plotting computations are floating illustrations.

The functions near the end of `run_paper_figures.py` contain one block per figure.
Colors, axis ranges, labels and display factors are explicit. Tests check spline
interpolation/norms, targets, reference evolutions, scales, legends and clipping.
They do not certify non-reachability.
