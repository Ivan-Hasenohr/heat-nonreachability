# Paper figures

Generate the seven figures from one completed campaign:

```console
python -m pip install -r figures/requirements.txt
python figures/run_paper_figures.py --project "results/campaigns/<timestamp>/code"
```

| Figure | Result file in the campaign's code/results/ | Figure filename |
|---|---|---|
| 1 | heat_intro.mat | heat_intro.pdf |
| 2 | heat_positive.mat | heat_positive.pdf |
| 3 | heat_min_time.mat | heat_min_time.pdf |
| 4 | heat_projected.mat | heat_projected.pdf |
| 5 | coupled_projected/ball/fine.mat | coupled_ball.pdf |
| 6 | coupled_projected/partial/fine.mat | coupled_partial.pdf |
| 7 | coupled_projected/split/fine.mat | coupled_split.pdf |

PNG and SVG versions are generated alongside each PDF. Output defaults to
`PROJECT/results/figures/` and can be changed with `--output PATH`.
The program checks all seven result files before generating figures.

Targets are plotted using the midpoints of saved coefficient enclosures.
Witnesses are natural splines reconstructed from the saved nodal values.
Both components of a coupled witness use one display factor. The partial
target's free component is explicitly labelled. Ball radii are product-L2
radii. The reference endpoints and all plotted values are floating
illustrations. `figure_metadata.json` records source hashes and display factors.
