"""IsoDec sequence matching in the wx GUI."""

import os
import subprocess
import sys
import unittest
from pathlib import Path


UNIDEC_ROOT = Path(__file__).resolve().parents[1]
@unittest.skipUnless(
    sys.platform.startswith("win") or sys.platform == "darwin"
    or bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")),
    "wxPython requires a graphical display",
)
class TestIsoDecFragmentWorkflow(unittest.TestCase):
    def test_sequence_panel_matches_current_peaks(self):
        script = """
import os
import sys
from types import SimpleNamespace
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.argv = ['isodec']
import isogen
import pandas as pd
from matplotlib.figure import Figure
from unidec.IsoDecGUI import IsoDecPres
import isodec
from isodec.match import MatchedCollection
from isodec.fragment_view import plot_fragment_matches

local_isodec = Path.cwd().parent / 'IsoDec' / 'isodec'
if (local_isodec / '__init__.py').is_file():
    assert Path(isodec.__file__).resolve().parent == local_isodec.resolve()

table = pd.DataFrame(index=range(1, 8), columns=['c_match', "z'_match"])
table.loc[4] = [1.0, 1.0]
ax = Figure().subplots()
plot_fragment_matches(ax, 'PEPTIDEK', SimpleNamespace(
    fragment_matches=table, sequence_coverage=1 / 7, fragment_match_percent=100),
    residues_per_line=4)
assert tuple(ax.lines[0].get_xdata()) == (3.3, 4)
assert tuple(ax.lines[2].get_xdata()) == (0, 0.7)

app = IsoDecPres()
try:
    view = app.view
    controls = view.controls
    assert view.sizerplot.GetItemSpan(view.fragment_panel) == (1, 2)
    assert controls.ctlfragmentation.GetValue() == 'ETD'
    assert controls.ctlfragmentppm.GetValue() == '5'
    assert controls.ctlmultiplemonoisotopics.GetValue()
    with patch('wx.MessageBox') as message:
        app.on_match_sequence()
    assert 'Run IsoDec' in message.call_args.args[0]

    mass = isogen.calc_pep_fragments('PEPTIDE', fragmentation_type='HCD')['b2']
    app.isodeceng.pks = MatchedCollection().add_peaks([SimpleNamespace(monoiso=mass)])
    app.match_source_ready = True
    controls.ctlsequence.SetValue('PEPTIDE')
    controls.ctlfragmentation.SetValue('HCD')
    controls.ctlfragmentppm.SetValue('5')
    app.on_match_sequence()
    assert app.isodeceng.pks.peaks[0].sequence_match == 'b2'
    assert len(view.fragment_ax.lines) == 2
    assert '16.7%' in view.GetStatusBar().GetStatusText(5)
    peak = app.isodeceng.pks.peaks[0]
    peak.monoiso = mass + 10
    peak.monoisos = [mass + 10, mass]
    controls.ctlmultiplemonoisotopics.SetValue(False)
    app.on_match_sequence()
    assert peak.sequence_match is None
    assert peak.monoisos == [mass + 10, mass]
    controls.ctlmultiplemonoisotopics.SetValue(True)
    app.on_match_sequence()
    assert peak.sequence_match == 'b2'
    controls.ctlavgpeakmasses.SetValue(True)
    with patch('unidec.IsoDecGUI.match_fragments', wraps=isodec.match_fragments) as matcher:
        app.on_match_sequence()
    assert matcher.call_args.kwargs['monoisotopic'] is False
    assert matcher.call_args.kwargs['match_multiple_monoisotopics'] is True
    controls.ctlavgpeakmasses.SetValue(False)
    app.on_match_sequence()
    with tempfile.TemporaryDirectory() as directory:
        flags, files = view.save_all_figures('pdf', header=directory + '/sequence')
        assert flags == [3]
        assert files[0][1].endswith('sequence_Figure3.pdf')
        assert os.path.isfile(files[0][1])

    controls.ctlfragmentppm.SetValue('-1')
    with patch('wx.MessageBox') as message:
        app.on_match_sequence()
    assert 'non-negative' in message.call_args.args[0]
    assert len(view.fragment_ax.lines) == 2

    with patch.object(app, 'on_dataprep_button'), patch.object(app, 'on_unidec_button'), patch.object(
        app, 'on_pick_peaks'
    ), patch.object(app, 'on_plot_peaks'), patch.object(app, 'on_plot_dists'), patch.object(
        app, 'on_match_sequence'
    ) as match_step:
        app.on_run_all()
        assert match_step.call_count == 1
        controls.ctlsequence.SetValue('')
        app.on_run_all()
        assert match_step.call_count == 1
    controls.ctlsequence.SetValue('PEPTIDE')

    with tempfile.TemporaryDirectory() as directory:
        old_dir, new_dir, empty_dir = (Path(directory) / name for name in ('old', 'new', 'empty'))
        for folder in (old_dir, new_dir, empty_dir):
            folder.mkdir()
        app.sequence_path = old_dir / 'seq.fasta'
        controls.ctlsequence.SetValue('S[Acetylation]HHS')
        (new_dir / 'seq.fasta').write_text('>previous sequence\\nS[Acetylation]\\nHHS\\n', encoding='utf-8')

        def open_into(folder):
            with patch.object(app, 'export_config'), patch.object(
                app.eng, 'open_file', side_effect=lambda *args, **kwargs: setattr(
                    app.eng.config, 'udir', str(folder))
            ), patch.object(app, 'makeplot1'), patch.object(app, 'import_config'), patch.object(
                app, 'write_to_recent'
            ), patch.object(view.menu, 'update_recent'):
                app.on_open_file('sample.txt', directory=directory)

        open_into(new_dir)
        assert (old_dir / 'seq.fasta').read_text(encoding='utf-8') == '>IsoDec sequence\\nS[Acetylation]HHS\\n'
        assert controls.ctlsequence.GetValue() == 'S[Acetylation]HHS'
        open_into(empty_dir)
        assert controls.ctlsequence.GetValue() == ''

        class StopProcess(Exception):
            pass

        controls.ctlsequence.SetValue('PEPTIDE')
        with patch.object(view, 'export_gui_to_config'), patch.object(
            app, 'fix_parameters'
        ), patch.object(app, 'translate_config'), patch.object(
            app, 'export_config'
        ), patch.object(app.isodeceng, 'batch_process_spectrum', side_effect=StopProcess):
            try:
                app.on_unidec_button()
            except StopProcess:
                pass
        assert (empty_dir / 'seq.fasta').read_text(encoding='utf-8') == '>IsoDec sequence\\nPEPTIDE\\n'

    view.clear_all_plots()
    assert not view.fragment_ax.lines
finally:
    app.view.Destroy()
    app.wx_app.Yield()
    app.wx_app.Destroy()
"""
        environment = os.environ.copy()
        environment["PYTHONPATH"] = os.pathsep.join(
            (str(UNIDEC_ROOT), environment.get("PYTHONPATH", ""))
        )
        result = subprocess.run(
            [sys.executable, "-c", script], cwd=UNIDEC_ROOT, env=environment,
            capture_output=True, text=True, timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
