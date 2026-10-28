import unittest
from types import SimpleNamespace

from matplotlib.figure import Figure

from unidec.modules.plotting.ZoomCommon import ZoomCommon


class TestPlottingLabels(unittest.TestCase):
    def setUp(self):
        self.figure = Figure()
        self.axes = self.figure.add_subplot()
        self.text = self.axes.text(1, 1, "label")
        self.line = self.axes.vlines(1, 0, 1)
        self.zoom = ZoomCommon()
        self.zoom.parent = SimpleNamespace(labeltexts=[self.text], labellines=[self.line])

    def test_kill_labels_removes_attached_artists(self):
        self.zoom.kill_labels(repaint=False)

        self.assertIsNone(self.text.axes)
        self.assertIsNone(self.line.axes)
        self.assertEqual(self.zoom.parent.labeltexts, [])
        self.assertEqual(self.zoom.parent.labellines, [])

    def test_kill_labels_ignores_artists_detached_by_figure_clear(self):
        self.figure.clear()

        self.zoom.kill_labels(repaint=False)

        self.assertEqual(self.zoom.parent.labeltexts, [])
        self.assertEqual(self.zoom.parent.labellines, [])


if __name__ == "__main__":
    unittest.main()
