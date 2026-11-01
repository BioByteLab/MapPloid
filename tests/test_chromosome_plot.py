from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import unittest

from mapploid.reporting import chromosome_support_figures, chromosome_rows, draw_chromosome_matrix, evidence_figure


class ChromosomePlotTests(unittest.TestCase):
    def test_empty_evidence_does_not_create_figures(self):
        with TemporaryDirectory() as temporary:
            out = Path(temporary) / 'figures'
            self.assertEqual(chromosome_support_figures({}, out), ([], []))
            self.assertFalse(out.exists())

    def test_matrix_orientation_uncertain_cells_and_natural_sequence_order(self):
        from matplotlib.figure import Figure
        rows = [dict(chrom='chr10', best_hypothesis='NA', status='TOO_FEW_SITES'),
                dict(chrom='chr2', best_hypothesis='4x', status='LOW_DELTA_BIC'),
                dict(chrom='chr1', best_hypothesis='2x', status='PASS'),
                dict(chrom='chr3', best_hypothesis='3x', status='UNRESOLVED')]
        ordered = chromosome_rows(dict(chromosomes=rows))
        ax = Figure().subplots()
        draw_chromosome_matrix(ax, ordered)
        self.assertEqual([t.get_text() for t in ax.get_xticklabels()], ['chr1', 'chr2', 'chr3', 'chr10'])
        self.assertEqual([t.get_text() for t in ax.get_yticklabels()], ['2x', '3x', '4x'])
        supported, candidate, unresolved, few = ax.patches
        self.assertAlmostEqual(supported.get_y(), .055)
        self.assertAlmostEqual(candidate.get_y(), 2.055)
        self.assertEqual(candidate.get_hatch(), '////')
        self.assertEqual(candidate.get_facecolor()[-1], 0)
        self.assertEqual((unresolved.get_y(), unresolved.get_height()), (0, 3))
        self.assertEqual((few.get_y(), few.get_height()), (0, 3))
        self.assertEqual([t.get_text() for t in ax.texts], ['NA', 'NA'])

    def test_multiple_pages_export_every_sequence_without_truncation(self):
        rows = [dict(chrom=f'chr{i}', sites=300, best_hypothesis='3x', delta_bic=20, status='PASS')
                for i in range(1, 26)]
        with TemporaryDirectory() as temporary:
            images, files = chromosome_support_figures(dict(chromosomes=rows), Path(temporary))
            self.assertEqual(len(images), 2)
            self.assertEqual([p.name for p in files], ['chromosome_support_001.pdf', 'chromosome_support_001.png',
                                                       'chromosome_support_002.pdf', 'chromosome_support_002.png'])
            self.assertTrue(all(p.stat().st_size > 0 for p in files))
            self.assertIn('>chr24</text>', images[0])
            self.assertNotIn('>chr25</text>', images[0])
            self.assertIn('>chr25</text>', images[1])

    def test_long_sequence_labels_do_not_overlap_at_full_column_count(self):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        from mapploid.reporting import chromosome_matrix_height, chromosome_matrix_panel
        rows = [dict(chrom=f'WWWW_chromosome_{i:02d}_long_scaffold_name',
                     best_hypothesis='3x', status='PASS') for i in range(12)]
        fig = Figure(figsize=(6.2, chromosome_matrix_height(rows)), layout='constrained')
        canvas = FigureCanvasAgg(fig)
        chromosome_matrix_panel(fig, rows)
        canvas.draw()
        labels = fig.axes[0].get_xticklabels()
        bounds = [label.get_window_extent(canvas.get_renderer()) for label in labels]
        for left, right in zip(bounds, bounds[1:]):
            self.assertLess(left.x1, right.x0)
        for label, row in zip(labels, rows):
            self.assertEqual(label.get_text().replace('\n', ''), row['chrom'])

    def test_combined_plot_contains_distributions_and_aligned_matrix(self):
        from mapploid import reporting
        stage = dict(model=dict(fits={str(p): dict(bic=p, converged=True) for p in (2,3,4)},
                                frequency_profile=dict(edges=[.15,.25,.35,.5], observed=[.2,.5,.3],
                                  fitted_models={str(p): [.2,.5,.3] for p in (2,3,4)})),
                     chromosomes=[dict(chrom='chr1', best_hypothesis='3x', status='PASS')])
        original = reporting.save_report_figure
        figures = []
        def capture(fig, *args):
            result = original(fig, *args)
            figures.append(fig)
            return result
        with patch.object(reporting, 'save_report_figure', capture):
            svg, _ = evidence_figure(stage)
        self.assertEqual(svg.count('<svg'), 1)
        self.assertIn('Chromosome model support', svg)
        self.assertTrue(all(f'{p}x model' in svg for p in (2,3,4)))
        self.assertEqual(len(figures[0].axes), 4)
        positions = [ax.get_position() for ax in figures[0].axes]
        self.assertLess(positions[3].y1, positions[2].y0)
        for box in positions[1:]:
            self.assertAlmostEqual(box.x0, positions[0].x0, places=5)
            self.assertAlmostEqual(box.x1, positions[0].x1, places=5)


if __name__ == '__main__':
    unittest.main()
