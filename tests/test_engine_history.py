import unittest
from types import SimpleNamespace

import numpy as np

from unidec.modules.unidec_enginebase import UniDecEngine


class TestEngineHistory(unittest.TestCase):
    def test_auto_peak_width_rejects_zero_data(self):
        engine = UniDecEngine(silent=True)
        engine.data = SimpleNamespace(data2=np.array([[100.0, 0.0], [101.0, 0.0]]))
        with self.assertRaisesRegex(ValueError, "positive finite intensity"):
            engine.get_auto_peak_width()

    def test_undo_redo_boundaries_and_snapshot_independence(self):
        engine = UniDecEngine(silent=True)
        initial_masslb = engine.config.masslb
        engine.clear_history()

        engine.config.masslb = 1234.0
        engine.config.smashlist = [[100.0, 2.0]]
        engine.update_history()
        engine.config.masslb = 5678.0
        engine.config.smashlist.append([200.0, 3.0])
        engine.update_history()

        self.assertEqual(engine.config_count, 3)
        self.assertEqual(engine.config_history[1].smashlist, [[100.0, 2.0]])

        engine.undo()
        self.assertEqual(engine.config.masslb, 1234.0)
        self.assertEqual(engine.config.smashlist, [[100.0, 2.0]])
        engine.undo()
        self.assertEqual(engine.config.masslb, initial_masslb)
        engine.undo()
        self.assertEqual(engine.config_count, 1)
        self.assertEqual(engine.config.masslb, initial_masslb)

        engine.redo()
        self.assertEqual(engine.config.masslb, 1234.0)
        engine.redo()
        self.assertEqual(engine.config.masslb, 5678.0)
        self.assertEqual(engine.config.smashlist, [[100.0, 2.0], [200.0, 3.0]])
        engine.redo()
        self.assertEqual(engine.config_count, 3)

        engine.config.smashlist.append([300.0, 4.0])
        self.assertEqual(engine.config_history[-1].smashlist, [[100.0, 2.0], [200.0, 3.0]])

    def test_edit_after_undo_discards_redo_branch(self):
        engine = UniDecEngine(silent=True)
        engine.clear_history()

        for masslb, smashlist in (
            (1000.0, [[100.0, 1.0]]),
            (2000.0, [[200.0, 2.0]]),
            (3000.0, [[300.0, 3.0]]),
        ):
            engine.config.masslb = masslb
            engine.config.smashlist = smashlist
            engine.update_history()

        engine.undo()
        self.assertEqual(engine.config.masslb, 2000.0)
        engine.config.masslb = 2500.0
        engine.config.smashlist.append([250.0, 2.5])
        engine.update_history()

        self.assertEqual(engine.config_count, 4)
        self.assertEqual(len(engine.config_history), 4)
        self.assertEqual(engine.config_history[-1].masslb, 2500.0)
        engine.redo()
        self.assertEqual(engine.config.masslb, 2500.0)
        self.assertEqual(engine.config_history[1].smashlist, [[100.0, 1.0]])
        self.assertEqual(engine.config_history[2].smashlist, [[200.0, 2.0]])


if __name__ == "__main__":
    unittest.main()
