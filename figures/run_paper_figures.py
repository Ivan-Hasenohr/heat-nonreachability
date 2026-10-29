"""Recreate the paper figures in Python, using the original notebook's style.

Run: python run_paper_figures.py
Refresh from MATLAB results: python run_paper_figures.py --project PATH
The bundled NPZ snapshots make subsequent runs independent of MATLAB/INTLAB.
All calculations here are floating illustrations, never CAP verification.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.io import loadmat
from scipy.linalg import expm


HERE = Path(__file__).resolve().parent
CASE_NAMES = ('corollary18', 'corollary19', 'proposition20', 'proposition22', 'proposition24')

# The visual choices follow Figures 1-3 and cells 6, 14 and 18 of the notebook.
STYLE = {
    'font.family': 'DejaVu Sans', 'font.size': 12,
    'mathtext.fontset': 'dejavusans', 'axes.linewidth': 0.8,
    'lines.linewidth': 1.5, 'axes.grid': False,
    'legend.fontsize': 12, 'legend.frameon': True,
    'legend.facecolor': 'white', 'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8', 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'svg.fonttype': 'none', 'savefig.facecolor': 'white',
}
INITIAL = 'green'
TARGET = 'red'
WITNESS = 'blue'
REGULARIZED = 'olive'
REFERENCE = 'lightgreen'
CONTROL = 'orange'


@dataclass
class Case:
    name: str
    p: np.ndarray
    nx: int
    nt: int
    time: float
    kappa: np.ndarray
    coupling: np.ndarray
    intervals: np.ndarray
    control_bounds: np.ndarray
    initial: np.ndarray
    target: np.ndarray
    certified: bool
    source: str
    source_sha256: str

    @property
    def components(self):
        return self.p.size // (self.nx - 1)

    def spline(self):
        interior = self.p.reshape(self.components, self.nx - 1).T
        nodes = np.vstack((np.zeros(self.components), interior, np.zeros(self.components)))
        return CubicSpline(np.linspace(0, 1, self.nx + 1), nodes, bc_type='natural')


def interval_bounds(value):
    """Read the old-style INTLAB MAT object as floating endpoints only."""
    fields = value.dtype.names or ()
    if 'inf' not in fields:
        numeric = np.asarray(value, dtype=float)
        return numeric, numeric
    assert value.size == 1, 'Expected one INTLAB object containing its array.'
    if np.any(value['complex'].item()):
        raise ValueError('Only real plotting data are supported.')
    lower, upper = value['inf'].item(), value['sup'].item()
    if np.size(lower) == 0:
        middle, radius = value['mid'].item(), value['rad'].item()
        lower, upper = middle - radius, middle + radius
    return np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)


def midpoint(value):
    lower, upper = interval_bounds(value)
    return (lower + upper) / 2


def read_case(project, name):
    raw = name == 'proposition20_unregularized'
    path = project / ('data/witnesses' if raw else 'results') / f'{name}.mat'
    names = ['p', 'args', 'y_0', 'y_0_base', 'y_f_base', 'certified',
             'certification_performed', 'lambda_reg']
    saved = loadmat(path, variable_names=names)
    args = saved['args'].ravel()
    p = midpoint(saved['p']).ravel()
    nx = int(midpoint(args[10]).item())
    components = p.size // (nx - 1)
    assert components in (1, 2) and p.size == components * (nx - 1)
    omega_index = 12 if components == 1 else 14
    lower, upper = interval_bounds(args[omega_index])
    intervals = np.column_stack((lower.ravel(), upper.ravel()))
    if raw:
        assert not saved['certification_performed'].item()
        assert saved['lambda_reg'].item() == 0
    else:
        assert saved['certified'].item(), f'{name} is not a saved certified result.'
    initial = saved.get('y_0_base', saved.get('y_0', np.zeros((components, 1))))
    return Case(
        name=name, p=p, nx=nx, nt=int(midpoint(args[1]).item()),
        time=float(midpoint(args[0]).item()), kappa=midpoint(args[11]).ravel(),
        coupling=midpoint(args[12]) if components == 2 else np.zeros((1, 1)),
        intervals=intervals,
        control_bounds=midpoint(args[13]).ravel() if components == 2
        else np.array([0., midpoint(args[8]).item()]),
        initial=midpoint(initial).ravel(), target=midpoint(saved['y_f_base']).ravel(),
        certified=not raw, source=path.relative_to(project).as_posix(), source_sha256=sha256(path),
    )


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_snapshot(case, folder):
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder / f'{case.name}.npz', **case.__dict__)


def read_snapshot(folder, name):
    with np.load(folder / f'{name}.npz', allow_pickle=False) as data:
        values = {key: data[key] for key in data.files}
    for key in ('name', 'nx', 'nt', 'time', 'certified', 'source', 'source_sha256'):
        values[key] = values[key].item()
    return Case(**values)


def state_values(base, nx, components, x):
    """Preserve P1 targets; do not smooth the tent with a cubic interpolant."""
    if base.size == components * (nx - 1):
        values = base.reshape(components, nx - 1).T
        nodes = np.vstack((np.zeros(components), values, np.zeros(components)))
        mesh = np.linspace(0, 1, nx + 1)
        return np.column_stack([np.interp(x, mesh, nodes[:, j]) for j in range(components)])
    if base.size == 2 * nx * components:
        coefficients = base.reshape(components, 2 * nx).T
        nodes = np.vstack((np.zeros(components), coefficients[:nx-1], np.zeros(components)))
        mu = coefficients[nx-1:]
        j = np.minimum((nx * x).astype(int), nx - 1)
        t = (nx * x - j)[:, None]
        return ((2*t**3-3*t**2+1)*nodes[j] + (t**3-2*t**2+t)*mu[j]
                + (-2*t**3+3*t**2)*nodes[j+1] + (t**3-t**2)*mu[j+1])
    coefficients = base.reshape(components, -1).T
    modes = np.arange(1, len(coefficients) + 1)
    return np.sin(np.pi * x[:, None] * modes) @ coefficients


def natural_spline_norm(spline):
    """Integrate the square cell by cell, using dimensionless coordinates."""
    h = np.diff(spline.x)
    a = spline.c * h[None, :, None] ** np.arange(3, -1, -1)[:, None, None]
    norm_squared = sum(
        np.sum(h[:, None] * a[i] * a[j]) / (7-i-j)
        for i in range(4) for j in range(4)
    )
    assert norm_squared > 0
    return np.sqrt(norm_squared)


def maximal_control_response(case, x, modes=128):
    """Continuous heat response: exact piecewise-quadratic steady profile
    minus its floating sine-series evolution. This is not a CAP enclosure.
    """
    kappa = case.kappa.item()
    maximum = case.control_bounds[1]
    a, b = case.intervals.T
    slope = maximum * np.sum((1-a)**2 - (1-b)**2) / (2*kappa)
    steady = slope*x - maximum/(2*kappa) * np.sum(
        np.maximum(x[:, None]-a, 0)**2 - np.maximum(x[:, None]-b, 0)**2, axis=1)
    n = np.arange(1, modes+1)
    forcing = 2*maximum*np.sum(np.cos(np.pi*n[:, None]*a)-np.cos(np.pi*n[:, None]*b), axis=1)
    coefficients = forcing/(kappa*(np.pi*n)**3)
    coefficients *= np.exp(-kappa*(np.pi*n)**2*case.time)
    return steady - np.sin(np.pi*x[:, None]*n) @ coefficients


def uncontrolled_final_state(case, x):
    coefficients = case.initial.reshape(case.components, -1).T
    evolved = np.array([
        expm((case.coupling-(np.pi*n)**2*np.diag(case.kappa))*case.time) @ row
        for n, row in enumerate(coefficients, start=1)
    ])
    return np.sin(np.pi*x[:, None]*np.arange(1, len(evolved)+1)) @ evolved


def control_region(ax, intervals, bottom, top):
    patches = []
    for j, (a, b) in enumerate(intervals):
        patches.append(ax.fill_between([a, b], bottom, top, facecolor=CONTROL,
                                       alpha=0.5, linewidth=0,
                                       label=r'$\omega$' if j == 0 else '_nolegend_', zorder=0))
    return patches[0]


def paper_axes(ax, limits, y_step):
    ax.set_xlim(0, 1)
    ax.set_ylim(*limits)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_locator(MultipleLocator(y_step))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f' if y_step < .1 else '%.1f'))
    ax.tick_params(direction='out', length=3.5, width=.8)


def scalar_figure(case, x, *, raw=None):
    """Figures 1 and 2: same short labels, colors and legend placement."""
    fig, ax = plt.subplots(figsize=(8, 6.6))
    limits = (-.02, .048) if raw is not None else ((-.044, .048) if case.name == 'corollary19' else (-.032, .048))
    paper_axes(ax, limits, .01)
    region = control_region(ax, case.intervals, 0, limits[1])
    initial = ax.plot(x, state_values(case.initial, case.nx, 1, x), color=INITIAL, label=r'$y_0$')[0]
    target = state_values(case.target, case.nx, 1, x)[:, 0]
    target_line = ax.plot(x, target, color=TARGET, label=r'$y_f$')[0]
    spline = case.spline()
    # Same continuous display norm for both tent witnesses, not equal maxima.
    display_norm = .025/np.sqrt(3) if raw is not None else .29/np.sqrt(185)
    factor = display_norm/natural_spline_norm(spline)
    witness = ax.plot(x, factor*spline(x)[:, 0], color=REGULARIZED if raw else WITNESS,
                      label=r'$p_f^{\mathrm{reg}}$' if raw else r'$p_f$')[0]
    reference = ax.plot(x, maximal_control_response(case, x), color='black', label=r'$L_T\overline{u}$')[0]
    metadata = {'witness_display_factor': factor, 'display_L2_norm': display_norm}
    if raw is None:
        handles = [initial, region, target_line, witness, reference]
        ax.legend(handles=handles, loc='lower center', ncol=5,
                  handlelength=1.8, columnspacing=1.35, handletextpad=.55)
    else:
        raw_spline = raw.spline()
        raw_factor = display_norm/natural_spline_norm(raw_spline)
        raw_line = ax.plot(x, raw_factor*raw_spline(x)[:, 0], color=WITNESS, label=r'$p_f$')[0]
        # Matplotlib fills legend columns vertically, as in the original figure.
        handles = [initial, region, target_line, raw_line, witness, reference]
        ax.legend(handles=handles, loc='lower right', ncol=3, handlelength=1.8,
                  columnspacing=1.3, handletextpad=.55)
        metadata['raw_witness_display_factor'] = raw_factor
        metadata['raw_snapshot_contains_certificate'] = False
    fig.subplots_adjust(left=.10, right=.985, bottom=.07, top=.985)
    return fig, metadata


def coupled_ball_figure(case, x):
    """Figure 3: common scale and the notebook's enlarged uncontrolled state."""
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(15, 6))
    initial = state_values(case.initial, case.nx, 2, x)
    target = state_values(case.target, case.nx, 2, x)
    witness = case.spline()(x)
    evolved = uncontrolled_final_state(case, x)
    factor = 1/np.max(np.abs(witness))
    evolution_factor = .72/np.max(np.abs(evolved))
    for ax in axes:
        paper_axes(ax, (-1.2, 1.4), .5)
    left, right = axes
    left.plot(x, initial[:, 0], color=INITIAL, label=r'$y_{01}$')
    left.plot(x, target[:, 0], color=TARGET, label=r'$y_{f1}$')
    left.plot(x, factor*witness[:, 0], color=WITNESS, label=r'$p_1$')
    left.plot(x, evolution_factor*evolved[:, 0], color=REFERENCE, label=r'$(S_Ty_0)_1$')
    control_region(right, case.intervals, -1, 1.4)
    right.plot(x, factor*witness[:, 1], color=WITNESS, label=r'$p_2$')
    right.plot(x, evolution_factor*evolved[:, 1], color=REFERENCE, label=r'$(S_Ty_0)_2$')
    assert np.all(initial[:, 1] == target[:, 1])
    right.plot(x, target[:, 1], color=TARGET, label=r'$y_{02}=y_{f2}$')
    for ax in axes:
        ax.legend(loc='upper center', ncol=4, handlelength=1.8, columnspacing=1.1, handletextpad=.5)
    fig.subplots_adjust(left=.055, right=.99, bottom=.10, top=.98, wspace=.20)
    return fig, {'witness_display_factor': factor, 'uncontrolled_state_display_factor': evolution_factor,
                 'target_ball_radius': .01, 'target_ball_is_not_a_pointwise_band': True}


def coupled_partial_figure(case, x):
    """Same two-panel style; explicitly distinguish a free target from zero."""
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(15, 6))
    target = state_values(case.target, case.nx, 2, x)
    witness = case.spline()(x)
    assert np.all(witness[:, 1] == 0)
    factor = .028/np.max(np.abs(witness))
    for ax in axes:
        paper_axes(ax, (-.035, .035), .01)
    left, right = axes
    left.plot(x, np.zeros_like(x), color=INITIAL, label=r'$y_{01}$')
    left.plot(x, target[:, 0], color=TARGET, label=r'$y_1$')
    left.plot(x, factor*witness[:, 0], color=WITNESS, label=r'$p_1$')
    control_region(right, case.intervals, -.035, .035)
    right.plot(x, witness[:, 1], color=WITNESS, label=r'$y_{02}=p_2=0$')
    right.text(.70, .75, r'$\mathcal{Y}_{f,2}=L^2(0,1)$', transform=right.transAxes,
               ha='center', va='center', fontsize=15)
    for ax in axes:
        ax.legend(loc='upper center', ncol=3, handlelength=1.8, columnspacing=1.2)
    fig.subplots_adjust(left=.065, right=.99, bottom=.10, top=.98, wspace=.20)
    return fig, {'witness_display_factor': factor, 'second_target_component': 'free'}


def export(fig, name, metadata, output):
    output.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(output / f'{name}.{extension}', dpi=220, bbox_inches='tight', pad_inches=.035)
    plt.close(fig)
    return {'case': name, **metadata}


def write_captions(records, output):
    """Keep article captions consistent with every regenerated display scale."""
    captions = []
    for record in records:
        name = record['case']
        factor = record['witness_display_factor']
        if name.startswith('corollary'):
            number = name.removeprefix('corollary')
            text = (rf'Target, control domain and dual certificate for Corollary~{number}. '
                    rf'The displayed certificate is ${factor:.6g}\,p_f$.')
        elif name == 'proposition20':
            raw_factor = record['raw_witness_display_factor']
            text = (r'Target and dual witnesses for Proposition~20. '
                    rf'The displayed curves are ${raw_factor:.6g}\,p_f$ and '
                    rf'${factor:.6g}\,p_f^{{\mathrm{{reg}}}}$, with the same displayed $L^2$ norm. '
                    r'The raw curve comes from a search snapshot; its fine-mesh verification is separate.')
        elif name == 'proposition22':
            state_factor = record['uncontrolled_state_display_factor']
            text = (r'Initial state, target centre, control domain and dual certificate for Proposition~22. '
                    rf'The displayed witness and uncontrolled final state are ${factor:.6g}\,p_f$ '
                    rf'and ${state_factor:.6g}\,S_Ty_0$, respectively. '
                    r'Each vector uses one factor for both components. '
                    r'The target is an $L^2$ ball, not a pointwise band.')
        else:
            text = (r'Prescribed first target component and dual certificate for Proposition~24. '
                    rf'The displayed witness is ${factor:.6g}\,p_f$. '
                    r'The second target component is free in $L^2(0,1)$, whereas $p_2=0$.')
        captions.append('\n'.join((r'\begin{figure}[htbp]', r'  \centering',
            rf'  \includegraphics[width=\linewidth]{{figures/{name}.pdf}}',
            rf'  \caption{{{text}}}', rf'  \label{{fig:{name}}}', r'\end{figure}')))
    (output/'captions.tex').write_text('\n\n'.join(captions)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, help='Read current MAT results from this paper-cap folder.')
    parser.add_argument('--data-dir', type=Path, default=HERE/'data', help='Portable floating snapshots.')
    parser.add_argument('--output', type=Path, default=HERE.parent/'results'/'figures')
    options = parser.parse_args()
    names = (*CASE_NAMES, 'proposition20_unregularized')
    use_snapshots = options.project is None and all((options.data_dir/f'{n}.npz').is_file() for n in names)
    if use_snapshots:
        cases = {n: read_snapshot(options.data_dir, n) for n in names}
    else:
        if options.project is None:
            parser.error('No bundled snapshots found. Supply --project PATH to paper-cap.')
        cases = {n: read_case(options.project, n) for n in names}
        for case in cases.values():
            save_snapshot(case, options.data_dir)
    x = np.unique(np.r_[np.linspace(0, 1, 4001), .2, .4, .5, .8])
    records = []
    with plt.rc_context(STYLE):
        # Corollary 18: y_0, y_f, omega, p_f and L_T u_bar.
        fig, meta = scalar_figure(cases['corollary18'], x)
        records.append(export(fig, 'corollary18', meta, options.output))

        # Corollary 19: its own saved witness and T = 0.13 response.
        fig, meta = scalar_figure(cases['corollary19'], x)
        records.append(export(fig, 'corollary19', meta, options.output))

        # Proposition 20: unregularized and regularized witnesses on one scale.
        fig, meta = scalar_figure(cases['proposition20'], x, raw=cases['proposition20_unregularized'])
        records.append(export(fig, 'proposition20', meta, options.output))

        # Proposition 22: same two-panel presentation as Figure 3.
        fig, meta = coupled_ball_figure(cases['proposition22'], x)
        records.append(export(fig, 'proposition22', meta, options.output))

        # Proposition 24: prescribed first target and free second component.
        fig, meta = coupled_partial_figure(cases['proposition24'], x)
        records.append(export(fig, 'proposition24', meta, options.output))

    provenance = {
        'figures': records, 'sources': {n: {'path': c.source, 'sha256': c.source_sha256,
            'Nx': c.nx, 'Nt': c.nt, 'T': c.time} for n, c in cases.items()},
        'notes': ['Floating illustrations only; no certification is run.',
                  'Witness and uncontrolled-state display factors do not alter saved proof data.',
                  'The raw tent snapshot contains search data; any fine-mesh verification is a separate computation.',
                  'Short labels follow the original paper; display scaling is documented separately.'],
    }
    (options.output/'figure_metadata.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    write_captions(records, options.output)
    print(f'Five figures exported as PDF, SVG and PNG to {options.output}')
    for record in records:
        print(record)


if __name__ == '__main__':
    main()
