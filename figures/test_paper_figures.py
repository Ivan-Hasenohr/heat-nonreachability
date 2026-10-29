"""Numerical and figure-contract checks; no claim of CAP validation."""
import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import run_paper_figures as paper


class PaperFigures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        names = (*paper.CASE_NAMES, 'proposition20_unregularized')
        cls.cases = {name: paper.read_snapshot(paper.HERE/'data', name) for name in names}
        cls.x = np.linspace(0, 1, 4001)

    def tearDown(self):
        plt.close('all')

    def test_parameters_match_the_five_cases(self):
        for name, time in zip(paper.CASE_NAMES, (1, .13, 1, .32, 1)):
            case = self.cases[name]
            self.assertAlmostEqual(case.time, time)
            self.assertTrue(case.certified)
            expected = [[.2, .4], [.8, 1]] if case.components == 1 else [[0, .5]]
            np.testing.assert_allclose(case.intervals, expected, atol=1e-15)
        self.assertFalse(self.cases['proposition20_unregularized'].certified)

    def test_spline_interpolates_the_frozen_nodes_and_zero_boundaries(self):
        for case in self.cases.values():
            spline = case.spline()
            mesh = np.linspace(0, 1, case.nx + 1)
            expected = np.vstack((np.zeros(case.components),
                                  case.p.reshape(case.components, case.nx-1).T,
                                  np.zeros(case.components)))
            np.testing.assert_allclose(spline(mesh), expected, atol=1e-12)
            np.testing.assert_allclose(spline([0, 1], 2), 0, atol=2e-8)

    def test_matches_the_cap_hermite_derivative_system(self):
        case = self.cases['proposition20_unregularized']
        nx = case.nx
        nodes = np.r_[0, case.p, 0]
        band = np.ones((3, nx+1))
        band[1] = 4
        band[1, [0, -1]] = 2
        rhs = 3*np.r_[nodes[1]-nodes[0], nodes[2:]-nodes[:-2], nodes[-1]-nodes[-2]]
        mu = solve_banded((1, 1), band, rhs)
        np.testing.assert_allclose(case.spline()(np.linspace(0, 1, nx+1), 1)[:, 0]/nx, mu, atol=1e-13)

    def test_exact_spline_norm_against_independent_gauss_quadrature(self):
        u, weights = leggauss(4)
        for case in self.cases.values():
            spline = case.spline()
            mesh = spline.x
            h = np.diff(mesh)
            samples = mesh[:-1, None]+h[:, None]*(u+1)/2
            squares = np.sum(spline(samples)**2, axis=-1)
            expected = np.sqrt(np.sum(squares*weights*h[:, None]/2))
            self.assertAlmostEqual(paper.natural_spline_norm(spline), expected, places=12)

    def test_tent_is_not_smoothed(self):
        case = self.cases['proposition20']
        values = paper.state_values(case.target, case.nx, 1, self.x)[:, 0]
        np.testing.assert_allclose(values, np.minimum(self.x/32, (1-self.x)/8), atol=1e-15)
        self.assertEqual(case.nx % 5, 0)

    def test_maximal_response_against_a_different_fourier_formula(self):
        x = np.unique(np.r_[np.linspace(0, 1, 129), .2, .4, .8])
        n = np.arange(1, 4097)
        for name in ('corollary18', 'corollary19', 'proposition20'):
            case = self.cases[name]
            a, b = case.intervals.T
            forcing = 2*np.sum(np.cos(np.pi*n[:, None]*a)-np.cos(np.pi*n[:, None]*b), axis=1)/(np.pi*n)
            decay = case.kappa.item()*(np.pi*n)**2
            reference = np.sin(np.pi*x[:, None]*n) @ (forcing*(-np.expm1(-decay*case.time))/decay)
            actual = paper.maximal_control_response(case, x)
            np.testing.assert_allclose(actual, reference, atol=1e-9)
            self.assertGreaterEqual(np.min(actual), -1e-15)
            np.testing.assert_allclose(actual[[0, -1]], 0, atol=1e-15)

    def test_coupled_evolution_against_closed_form_two_by_two_exponential(self):
        case = self.cases['proposition22']
        matrix = case.coupling-np.pi**2*np.diag(case.kappa)
        center = np.trace(matrix)/2
        shifted = matrix-center*np.eye(2)
        delta = np.sqrt(shifted[0, 0]**2+shifted[0, 1]*shifted[1, 0])
        exponential = np.exp(center*case.time)*(np.cosh(delta*case.time)*np.eye(2)
                       + np.sinh(delta*case.time)/delta*shifted)
        expected = np.sin(np.pi*self.x[:, None])*(exponential @ np.array([1, 0]))
        np.testing.assert_allclose(paper.uncontrolled_final_state(case, self.x), expected, atol=1e-15)

    def test_tent_witnesses_have_one_common_display_norm(self):
        regularized = self.cases['proposition20']
        raw = self.cases['proposition20_unregularized']
        before = raw.p.copy()
        _, metadata = paper.scalar_figure(regularized, self.x, raw=raw)
        expected = .025/np.sqrt(3)
        self.assertAlmostEqual(metadata['witness_display_factor']*paper.natural_spline_norm(regularized.spline()), expected)
        self.assertAlmostEqual(metadata['raw_witness_display_factor']*paper.natural_spline_norm(raw.spline()), expected)
        self.assertFalse(metadata['raw_snapshot_contains_certificate'])
        np.testing.assert_array_equal(raw.p, before)

    def test_colors_and_short_notation_follow_the_paper(self):
        fig, _ = paper.scalar_figure(self.cases['corollary18'], self.x)
        ax = fig.axes[0]
        colors = {line.get_label(): to_rgba(line.get_color()) for line in ax.lines}
        self.assertEqual(colors[r'$y_0$'], to_rgba('green'))
        self.assertEqual(colors[r'$y_f$'], to_rgba('red'))
        self.assertEqual(colors[r'$p_f$'], to_rgba('blue'))
        self.assertEqual(colors[r'$L_T\overline{u}$'], to_rgba('black'))
        self.assertEqual([t.get_text() for t in ax.get_legend().texts],
                         [r'$y_0$', r'$\omega$', r'$y_f$', r'$p_f$', r'$L_T\overline{u}$'])

    def test_partial_target_is_free_and_not_a_zero_target_curve(self):
        fig, metadata = paper.coupled_partial_figure(self.cases['proposition24'], self.x)
        self.assertEqual(metadata['second_target_component'], 'free')
        self.assertTrue(any(r'L^2' in text.get_text() for text in fig.axes[1].texts))
        self.assertFalse(any(line.get_label() == r'$y_{f2}$' for line in fig.axes[1].lines))
        np.testing.assert_array_equal(self.cases['proposition24'].spline()(self.x)[:, 1], 0)

    def test_common_component_scaling_and_correct_control_component(self):
        fig, metadata = paper.coupled_ball_figure(self.cases['proposition22'], self.x)
        self.assertEqual(fig.axes[0].get_ylim(), fig.axes[1].get_ylim())
        self.assertEqual(len(fig.axes[0].collections), 0)
        self.assertEqual(len(fig.axes[1].collections), 1)
        self.assertGreater(metadata['uncontrolled_state_display_factor'], 1)
        self.assertTrue(metadata['target_ball_is_not_a_pointwise_band'])

    def test_no_clipped_curves_or_outside_legends_or_diagnostic_titles(self):
        figures = [paper.scalar_figure(self.cases[n], self.x)[0] for n in ('corollary18', 'corollary19')]
        figures.append(paper.scalar_figure(self.cases['proposition20'], self.x,
                                          raw=self.cases['proposition20_unregularized'])[0])
        figures.extend((paper.coupled_ball_figure(self.cases['proposition22'], self.x)[0],
                        paper.coupled_partial_figure(self.cases['proposition24'], self.x)[0]))
        for fig in figures:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            for ax in fig.axes:
                self.assertEqual(ax.get_title(), '')
                lo, hi = ax.get_ylim()
                for line in ax.lines:
                    self.assertGreaterEqual(np.min(line.get_ydata()), lo)
                    self.assertLessEqual(np.max(line.get_ydata()), hi)
                axes_box = ax.get_window_extent(renderer)
                legend_box = ax.get_legend().get_window_extent(renderer)
                self.assertTrue(axes_box.contains(legend_box.x0, legend_box.y0))
                self.assertTrue(axes_box.contains(legend_box.x1, legend_box.y1))


if __name__ == '__main__':
    unittest.main(verbosity=2)
