"""Independent plotting checks using synthetic inputs, without campaign results."""
import tempfile
import unittest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.io import savemat
from scipy.linalg import expm, solve_banded

import run_paper_figures as paper


def fixture(name, components=1):
    nx = 20
    x = np.arange(1, nx) / nx
    p = np.sin(np.pi*x)
    target = np.array([0.035])
    kappa, coupling = np.array([1.0]), np.zeros((1, 1))
    intervals = np.array([[0.2, 0.4], [0.8, 1.0]])
    bounds, kind, radius = np.array([0.0, 1.0]), 2, 0.0
    if components == 2:
        p = np.r_[p, np.zeros_like(p) if name == 'coupled_partial' else np.sin(2*np.pi*x)]
        target = np.array([0.001, -0.0005, 0.002, 0.001])
        kappa, coupling = np.array([1.0, 1.5]), np.array([[-1, 2], [-2, -1]])
        intervals, bounds = np.array([[0.0, 0.5]]), np.array([-1.0, 2.0])
        kind = 1 if name == 'coupled_partial' else 2
        radius = 0.0002 if name == 'coupled_ball' else 0.0
    return paper.Case(name, p, nx, 100, 0.13, kappa, coupling, intervals, bounds,
                      np.zeros(components), target, False, 'synthetic fixture', '', kind, radius)


class PaperFigures(unittest.TestCase):
    def setUp(self):
        self.cases = [fixture('heat_positive'), fixture('coupled_ball', 2),
                      fixture('coupled_partial', 2), fixture('coupled_split', 2)]
        self.x = np.linspace(0, 1, 2001)

    def tearDown(self):
        plt.close('all')

    def test_spline_interpolates_nodes_and_boundaries(self):
        for case in self.cases:
            expected = np.pad(case.p.reshape(case.components, -1).T, ((1, 1), (0, 0)))
            spline = case.spline()
            np.testing.assert_allclose(spline(np.linspace(0, 1, case.nx+1)), expected, atol=1e-14)
            np.testing.assert_allclose(spline([0, 1], 2), 0, atol=1e-11)

    def test_spline_derivatives_against_hermite_system(self):
        case = self.cases[0]
        nodes = np.r_[0, case.p, 0]
        band = np.ones((3, case.nx+1))
        band[1] = 4
        band[1, [0, -1]] = 2
        rhs = 3*np.r_[nodes[1]-nodes[0], nodes[2:]-nodes[:-2], nodes[-1]-nodes[-2]]
        mu = solve_banded((1, 1), band, rhs)
        actual = case.spline()(np.linspace(0, 1, case.nx+1), 1)[:, 0]/case.nx
        np.testing.assert_allclose(actual, mu, atol=1e-13)

    def test_spline_norm_against_gauss_quadrature(self):
        points, weights = leggauss(4)
        for case in self.cases:
            spline = case.spline()
            h = np.diff(spline.x)
            samples = spline.x[:-1, None] + h[:, None]*(points+1)/2
            expected = np.sqrt(np.sum(np.sum(spline(samples)**2, axis=-1)*weights*h[:, None]/2))
            self.assertAlmostEqual(paper.natural_spline_norm(spline), expected, places=12)

    def test_scalar_reference_against_fourier_formula(self):
        case = self.cases[0]
        x = np.linspace(0, 1, 129)
        k = np.arange(1, 4097)
        a, b = case.intervals.T
        forcing = 2*np.sum(np.cos(np.pi*k[:, None]*a)-np.cos(np.pi*k[:, None]*b), axis=1)/(np.pi*k)
        decay = (np.pi*k)**2
        expected = np.sin(np.pi*x[:, None]*k) @ (forcing*(-np.expm1(-decay*case.time))/decay)
        np.testing.assert_allclose(paper.maximal_control_response(case, x), expected, atol=1e-9)

    def test_coupled_reference_against_matrix_exponential(self):
        profile = [(0, 0.25, 1), (0.25, 0.5, -0.5)]
        actual = paper.endpoint_coefficients(0.06, profile, 4)
        for j in range(4):
            k = j+1
            lam = (np.pi*k)**2
            generator = np.array([[-lam-1, 2], [-2, -1.5*lam-1]])
            forcing = sum(2*value*(np.cos(np.pi*k*a)-np.cos(np.pi*k*b))/(np.pi*k)
                          for a, b, value in profile)
            expected = np.linalg.solve(generator, (expm(0.06*generator)-np.eye(2)) @ [0, forcing])
            np.testing.assert_allclose(actual[j], expected, rtol=1e-12, atol=1e-16)

    def test_coupled_scale_and_control_component(self):
        case = self.cases[1]
        fig, metadata = paper.coupled_figure(case, self.x)
        self.assertEqual(fig.axes[0].get_ylim(), fig.axes[1].get_ylim())
        self.assertEqual(len(fig.axes[0].collections), 0)
        self.assertEqual(len(fig.axes[1].collections), 1)
        self.assertTrue(metadata['shared_component_scale'])
        self.assertEqual(metadata['target_ball_radius'], case.target_radius)

    def test_free_component_and_zero_witness(self):
        fig, metadata = paper.coupled_figure(self.cases[2], self.x)
        self.assertTrue(metadata['free_second_component'])
        self.assertTrue(any('L^2' in text.get_text() for text in fig.axes[1].texts))
        np.testing.assert_array_equal(self.cases[2].spline()(self.x)[:, 1], 0)

    def test_visible_curves_are_inside_axes(self):
        for case in self.cases:
            fig, _ = paper.scalar_figure(case, self.x) if case.components == 1 else paper.coupled_figure(case, self.x)
            for ax in fig.axes:
                lo, hi = ax.get_ylim()
                for line in ax.lines:
                    self.assertGreaterEqual(float(np.min(line.get_ydata())), lo)
                    self.assertLessEqual(float(np.max(line.get_ydata())), hi)

    def test_loader_rejects_nonnegative_or_uncertified_results(self):
        case = self.cases[0]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'results').mkdir()
            data = dict(p=case.p, Nx=case.nx, Nt=case.nt, T=case.time, kappa=case.kappa,
                        control_intervals=np.array([[0.2, 0.4]]), control_bounds=case.control_bounds,
                        y_f_base=case.target, certified=True, J=0.0)
            source = root/'results'/'heat_positive.mat'
            savemat(source, data)
            with self.assertRaises(ValueError):
                paper.read_case(root, 'heat_positive')
            data.update(J=-1.0, certified=False)
            savemat(source, data)
            with self.assertRaises(ValueError):
                paper.read_case(root, 'heat_positive')


if __name__ == '__main__':
    unittest.main(verbosity=2)
