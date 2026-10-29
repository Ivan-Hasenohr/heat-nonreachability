# Paper figures

The paper figures use the Python generator documented in
[../figures/README.md](../figures/README.md).

From the repository root:

```console
python -m pip install -r figures/requirements.txt
python figures/run_paper_figures.py
python figures/test_paper_figures.py
```

The default inputs are portable floating snapshots of the reference campaign.
To refresh them from a completed campaign, pass its `code` directory:

```console
python figures/run_paper_figures.py --project "results/campaigns/<timestamp>/code"
```

No optimization or certification is run by the plotting code. The former MATLAB
paper-plotting module has been removed; the two basic functions remaining in
`src/plotting` are unrelated optional helpers.
