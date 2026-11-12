"""Generate the seven paper figures from one completed certification campaign."""

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

# Shared style for the seven figures in the article.
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
    target_kind: int = 2
    target_radius: float = 0.0

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


# Stable physical names, ordered as Figures 1 to 7 in the manuscript.
RESULT_FILES = {
    'heat_intro': 'heat_intro.mat',
    'heat_positive': 'heat_positive.mat',
    'heat_min_time': 'heat_min_time.mat',
    'heat_projected': 'heat_projected.mat',
    'coupled_ball': 'coupled_projected/ball/fine.mat',
    'coupled_partial': 'coupled_projected/partial/fine.mat',
    'coupled_split': 'coupled_projected/split/fine.mat',
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_case(project, name):
    path = project / 'results' / RESULT_FILES[name]
    variables = ['p', 'Nx', 'Nt', 'T', 'kappa', 'reaction', 'control_intervals',
                 'control_bounds', 'y_0', 'y_0_base', 'y_f_base', 'target_kind',
                 'target_radius', 'certified', 'J']
    saved = loadmat(path, variable_names=variables)
    lower_j, upper_j = interval_bounds(saved['J'])
    if not saved['certified'].item() or not (
            np.isfinite(lower_j.item()) and np.isfinite(upper_j.item()) and upper_j.item() < 0):
        raise ValueError(f'{path} does not contain a finite negative certificate.')
    nx = int(saved['Nx'].item())
    p = saved['p'].ravel()
    components = p.size // (nx - 1)
    if components not in (1, 2) or p.size != components * (nx - 1):
        raise ValueError(f'{path} has inconsistent witness dimensions.')
    left, right = interval_bounds(saved['control_intervals'])
    initial = saved.get('y_0_base', saved.get('y_0', np.zeros((components, 1))))
    return Case(
        name=name, p=p, nx=nx, nt=int(saved['Nt'].item()), time=midpoint(saved['T']).item(),
        kappa=midpoint(saved['kappa']).ravel(),
        coupling=midpoint(saved.get('reaction', np.zeros((1, 1)))),
        intervals=np.column_stack((left.ravel(), right.ravel())),
        control_bounds=midpoint(saved['control_bounds']).ravel(),
        initial=midpoint(initial).ravel(), target=midpoint(saved['y_f_base']).ravel(),
        certified=True, source=path.relative_to(project).as_posix(), source_sha256=sha256(path),
        target_kind=int(saved.get('target_kind', np.array([[2]])).item()),
        target_radius=midpoint(saved.get('target_radius', np.array([[0.0]]))).item(),
    )


def scalar_figure(case, x):
    fig, ax = plt.subplots(figsize=(8, 6.6))
    spline = case.spline()
    display_norm = 0.29 / np.sqrt(185)
    factor = display_norm / natural_spline_norm(spline)
    witness = factor * spline(x)[:, 0]
    target = state_values(case.target, case.nx, 1, x)[:, 0]
    reference_level = 0.75 if case.name == 'heat_projected' else 1.0
    reference = reference_level * maximal_control_response(case, x)
    bottom = min(-0.032, float(witness.min()) - 0.005, float(target.min()) - 0.005)
    top = max(0.048, float(witness.max()) + 0.005, float(reference.max()) + 0.005)
    paper_axes(ax, (bottom, top), 0.01)
    region = control_region(ax, case.intervals, 0, top)
    initial = ax.plot(x, np.zeros_like(x), color=INITIAL, label=r'$y_0$')[0]
    target_line = ax.plot(x, target, color=TARGET, label=r'$y_f$')[0]
    witness_line = ax.plot(x, witness, color=WITNESS, label=r'$p_f$')[0]
    reference_line = ax.plot(x, reference, color='black', linestyle=(0, (2, 2)),
                             linewidth=0.9, label=r'$L_T\overline{u}$')[0]
    ax.legend(handles=[initial, region, target_line, witness_line, reference_line],
              loc='lower center', ncol=5, handlelength=1.8, columnspacing=1.35, handletextpad=0.55)
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.07, top=0.985)
    return fig, {'witness_display_factor': factor, 'display_L2_norm': display_norm,
                 'reference_control_level': reference_level}


def endpoint_coefficients(time, profile, n_modes=2048):
    """Floating reference endpoint for the coupled operator used in the article."""
    frequency = np.pi * np.arange(1, n_modes + 1)
    lam = frequency**2
    alpha, beta = lam + 1, 1.5 * lam + 1
    difference = lam / 4
    delta = np.sqrt(difference**2 - 4)
    tau = -(alpha + beta) / 2
    slow, fast = np.exp(time * (tau + delta)), np.exp(time * (tau - delta))
    e12 = (slow - fast) / delta
    e22 = ((1 - difference / delta) * slow + (1 + difference / delta) * fast) / 2
    forcing = sum(2 * level * (np.cos(frequency * a) - np.cos(frequency * b)) / frequency
                  for a, b, level in profile)
    return (forcing / (alpha * beta + 4))[:, None] * np.column_stack((
        -beta * e12 - 2 * (e22 - 1), 2 * e12 - alpha * (e22 - 1)))


def coupled_figure(case, x):
    from matplotlib.ticker import FuncFormatter

    if not (np.array_equal(case.kappa, [1, 1.5])
            and np.array_equal(case.coupling, [[-1, 2], [-2, -1]])):
        raise ValueError('The plotted coupled reference requires the paper model.')
    target = state_values(case.target, case.nx, 2, x)
    witness = case.spline()(x)
    factor = np.max(np.abs(target)) / np.max(np.abs(witness))
    displayed = factor * witness
    reference = None
    if case.name != 'coupled_partial':
        profile = [(0, 0.5, 1.8)] if case.name == 'coupled_ball' else [(0, 0.25, 1), (0.25, 0.5, -0.5)]
        coefficients = endpoint_coefficients(case.time, profile)
        reference = np.sin(np.pi * np.outer(x, np.arange(1, len(coefficients) + 1))) @ coefficients
    visible = [target, displayed] + ([] if reference is None else [reference])
    low = min(0, *(float(values.min()) for values in visible))
    high = max(0, *(float(values.max()) for values in visible))
    span = high - low
    limits = (low - 0.12 * span, high + 0.29 * span)
    wanted = (limits[1] - limits[0]) / 7
    choices = np.array([1, 2, 5, 10]) * 10**np.floor(np.log10(wanted))
    step = float(choices[np.argmin(np.abs(np.log(choices / wanted)))])
    precision = max(0, int(-np.floor(np.log10(step))))
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(15, 6))
    for component, ax in enumerate(axes):
        paper_axes(ax, limits, step)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, position: f'{value:.{precision}f}'))
        ax.set_title('First component' if component == 0 else 'Second component', pad=12)
        if component == 1:
            control_region(ax, case.intervals, limits[0], limits[1])
        ax.plot(x, np.zeros_like(x), color=INITIAL, label=rf'$y_{{0,{component+1}}}=0$')
        if case.target_kind == 1 and component == 1:
            ax.text(0.72, 0.72, r'$\mathcal{Y}_{f,2}=L^2(0,1)$', transform=ax.transAxes,
                    ha='center', fontsize=15)
            ax.text(0.72, 0.62, 'Free target component', transform=ax.transAxes, ha='center', fontsize=11)
        else:
            label = rf'$y_{{f,{component+1}}}$'
            if case.target_radius > 0:
                label = 'Center ' + label
            ax.plot(x, target[:, component], color=TARGET, label=label)
        ax.plot(x, displayed[:, component], color=WITNESS,
                linestyle='--' if np.all(displayed[:, component] == 0) else '-',
                label=rf'$({factor:.3g})\,p_{{f,{component+1}}}$')
        if reference is not None:
            ax.plot(x, reference[:, component], color='black', linestyle=(0, (2, 2)),
                    linewidth=0.9, label=r'$L_Tu_{\mathrm{ref}}$')
        ax.legend(loc='upper center', ncol=2, handlelength=1.8, columnspacing=1.2, handletextpad=0.5)
    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.10, top=0.89, wspace=0.16)
    return fig, {'witness_display_factor': factor, 'shared_component_scale': True,
                 'target_ball_radius': case.target_radius,
                 'free_second_component': case.target_kind == 1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True,
                        help='Completed campaign code directory or repository with standalone results.')
    parser.add_argument('--output', type=Path, help='Defaults to PROJECT/results/figures.')
    options = parser.parse_args()
    project = options.project.resolve()
    output = options.output or project / 'results/figures'
    # Load every input before writing any figure, so a missing result fails early.
    cases = [read_case(project, name) for name in RESULT_FILES]
    output.mkdir(parents=True, exist_ok=True)
    records = []
    x = np.unique(np.r_[np.linspace(0, 1, 6001), 0.2, 0.25, 0.4, 0.5, 0.8])
    with plt.rc_context(STYLE):
        for number, case in enumerate(cases, 1):
            figure, record = scalar_figure(case, x) if case.components == 1 else coupled_figure(case, x)
            for extension in ('pdf', 'svg', 'png'):
                figure.savefig(output / f'{case.name}.{extension}', dpi=220,
                               bbox_inches='tight', pad_inches=0.035)
            plt.close(figure)
            record.update(figure=number, case=case.name, source=case.source,
                          source_sha256=case.source_sha256, Nx=case.nx, Nt=case.nt, T=case.time,
                          witness_representation='natural cubic spline', floating_illustration=True)
            records.append(record)
    (output / 'figure_metadata.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
    print(f'Seven figures generated from {project} in {output}')


if __name__ == '__main__':
    main()
