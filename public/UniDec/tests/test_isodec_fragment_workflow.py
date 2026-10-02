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
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from unidec.IsoDecGUI import IsoDecPres
import wx
import isodec
from isodec.match import MatchedCollection
from isodec.isotope import calc_isotope_dist_dual
from isodec.fragment_view import plot_fragment_matches
from unidec.modules.peakstructure import Peaks

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

recent_path = os.path.join(tempfile.gettempdir(), 'isodec-recent-test.txt')
with patch.object(IsoDecPres, 'read_recent', return_value=[recent_path]):
    app = IsoDecPres()
try:
    view = app.view
    view.SetSize((984, 728))
    view.Show()
    app.wx_app.Yield()
    assert view.fragment_panel.IsShown()
    assert view.plot2.GetPosition().x + view.plot2.GetSize().width <= view.plotpanel.GetClientSize().width
    assert all(plot.canvas.GetSize() == plot.GetSize() for plot in (view.plot1, view.plot2))
    assert view.fragment_canvas.GetSize() == view.fragment_panel.GetSize()
    assert view.plotpanel.GetVirtualSize().height == view.plotpanel.GetClientSize().height
    assert view.controls.foldpanels.GetFoldPanel(6).IsExpanded()
    assert view.peakpanel.list_ctrl.GetSize().height >= view.peakpanel.GetClientSize().height - 4
    assert view.sizerplot.GetSize().height == view.plotpanel.GetVirtualSize().height
    if view.GetClientSize().width >= 1300:
        assert view.peakpanel.GetSize().width < 300
        assert view.plotpanel.GetSize().width >= 800
    for plot in (view.plot1, view.plot2):
        plot.centroid_plot(np.array([[500, 1.2e8], [600, 2.4e8], [700, 1.1e8]]),
                           xlabel='Mass', ylabel='Intensity')
        plot.canvas.draw()
        assert plot.subplot1.yaxis.get_tightbbox(plot.canvas.get_renderer()).x0 >= 0
        plot.clear_plot()
    recent_item = view.menu.menuOpenRecent.GetMenuItems()[0]
    assert recent_item.GetItemLabelText() == os.path.basename(recent_path)
    with patch.object(app, 'on_open_file') as open_recent:
        view.ProcessEvent(wx.CommandEvent(wx.EVT_MENU.typeId, recent_item.GetId()))
    open_recent.assert_called_once_with(os.path.basename(recent_path), os.path.dirname(recent_path))
    examples = [Path(path) for _, path, _ in view.menu.masterd2]
    ca_etd = next(path for path in examples if path.name == 'ca_etd.dat')
    assert (ca_etd.with_suffix('').with_name('ca_etd_unidecfiles') / 'seq.fasta').is_file()
    with patch.object(app, 'on_open_file') as open_example:
        view.menu.load_example_data(examples.index(ca_etd))
    open_example.assert_called_once_with(ca_etd.name, str(ca_etd.parent))
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
    assert view.fragment_panel.IsShown()
    assert len(view.fragment_ax.lines) == 2
    assert 'Fragments Matched' in view.fragment_ax.get_title(loc='left')
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
                assert app.eng.open_file.call_args.kwargs['simple_output'] is True

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
        ), patch.object(view.peakpanel, 'clear_list') as clear_peaks, patch.object(
            app.isodeceng, 'batch_process_spectrum', side_effect=StopProcess
        ):
            try:
                app.on_unidec_button()
            except StopProcess:
                pass
        clear_peaks.assert_called_once()
        assert (empty_dir / 'seq.fasta').read_text(encoding='utf-8') == '>IsoDec sequence\\nPEPTIDE\\n'

    view.clear_all_plots()
    assert not view.fragment_ax.lines

    assert controls.bruteforcebutton.GetLabel() == 'Brute Force Match'
    controls.ctlfragmentation.SetValue('HCD')
    controls.ctlfragmentppm.SetValue('5')
    controls.ctlsequence.SetValue('PEPTIDE')
    controls.ctlcentroided.SetValue(True)
    batch = isogen.calc_pep_fragment_isodists('PEPTIDE', fragmentation_type='HCD')
    index = batch.labels.index('b6')
    values = batch.intensities[index]
    positions = np.flatnonzero(values > values.max() * 0.01)
    spectrum = np.column_stack((batch.masses[index] / 2 + positions * 1.0033 / 2 + 1.007276467,
                                values[positions] * 100))
    app.eng.data.rawdata = spectrum
    app.eng.data.data2 = spectrum
    app.sequence_path = None
    controls.ctlminmz.SetValue(str(spectrum[0, 0] - 1))
    controls.ctlmaxmz.SetValue(str(spectrum[-1, 0] + 1))
    with patch('wx.MessageBox') as brute_message, patch.object(app, 'makeplot1'), patch.object(
        app, 'makeplot2'
    ):
        app.on_brute_force_match()
    assert brute_message.call_args is None
    assert any(peak.sequence_match == 'b6' and peak.z == 2 for peak in app.isodeceng.pks)
    assert app.isodeceng.pks.masses
    assert np.isclose(app.isodeceng.pks.fragment_matches.loc[6, 'b_match'], batch.masses[index])
    assert view.fragment_panel.IsShown()
    assert view.fragment_has_matches
    assert len(view.fragment_ax.lines) == 2
    assert view.fragment_ax.get_title(loc='left').startswith('Sequence Coverage:')
    assert 'Fragments Matched' not in view.fragment_ax.get_title(loc='left')
    assert view.peakpanel.list_ctrl.GetItemCount() == len(app.eng.pks.peaks)
    assert any('b6' in view.peakpanel.list_ctrl.GetItem(i, 4).GetText()
               for i in range(view.peakpanel.list_ctrl.GetItemCount()))
    assert 'Brute Force Match:' in view.GetStatusBar().GetStatusText(5)
    app.on_match_sequence()
    assert any('b6' in view.peakpanel.list_ctrl.GetItem(i, 4).GetText()
               for i in range(view.peakpanel.list_ctrl.GetItemCount()))
    controls.ctlavgpeakmasses.SetValue(True)
    with patch('wx.MessageBox') as message:
        app.on_brute_force_match()
    assert 'Turn off Average Mass' in message.call_args.args[0]

    peaks = Peaks()
    peaks.add_peaks(np.array([[600.12341, 10], [600.12349, 8]]))
    peaks.default_params()
    for index, peak in enumerate(peaks.peaks):
        peak.mztab = np.array([[500 + index, 10]])
        peak.stickdat = np.array([[500 + index, 5]])
        peak.avgmass = 601.12341 + index
    app.eng.pks = peaks
    app.eng.data.data2 = np.array([[500, 10], [501, 8]])
    app.eng.data.massdat = np.array([[600, 10], [601, 8]])
    view.peakpanel.add_data(peaks, collab1='Avg Mass')
    assert view.peakpanel.list_ctrl.GetItem(0, 1).GetText() == '601.123'
    assert view.peakpanel._item_mass(0) == peaks.peaks[0].mass
    view.peakpanel.add_data(peaks, show='zs')
    assert view.peakpanel.list_ctrl.GetColumnWidth(1) == 80
    assert view.peakpanel.list_ctrl.GetColumnWidth(3) == 40
    assert view.peakpanel.list_ctrl.GetColumnWidth(4) == 55
    assert [view.peakpanel.list_ctrl.GetItem(i, 1).GetText() for i in range(2)] == [
        '600.123', '600.123']
    app.plot_mass_peaks()
    app.plot_mz_peaks()
    app.on_plot_dists()
    assert sum(line.get_gid() == 'isodec_isotope' for line in view.plot1.subplot1.lines) == 2
    view.peakpanel.list_ctrl.Select(0)
    view.peakpanel.on_popup_two()
    assert [peak.ignore for peak in peaks.peaks] == [0, 1]
    isotopes = [line for line in view.plot1.subplot1.lines if line.get_gid() == 'isodec_isotope']
    assert len(isotopes) == 1
    assert 500 in isotopes[0].get_xdata()
    view.peakpanel.on_popup_three()
    assert sum(line.get_gid() == 'isodec_isotope' for line in view.plot1.subplot1.lines) == 2
    app.plot_mz_peaks()
    view.peakpanel.list_ctrl.Select(0)
    view.peakpanel.on_popup_two()
    assert not any(line.get_gid() == 'isodec_isotope' for line in view.plot1.subplot1.lines)

    with tempfile.TemporaryDirectory() as directory:
        app.eng.config.udir = directory
        app.eng.config.idconfig.write_msalign = 1
        app.eng.config.idconfig.write_tsv = 1
        with patch.object(app.isodeceng, 'export_peaks') as export:
            app.export_results()
        assert [call.kwargs['filename'] for call in export.call_args_list] == [
            str(Path(directory) / 'results'), str(Path(directory) / 'results.tsv')]

        spectrum_path = Path(directory) / 'sample.dat'
        spectrum_path.write_text('500 1\\n501 2\\n', encoding='utf-8')
        current_directory = os.getcwd()
        with patch.object(app.isodeceng, 'process_file'), patch.object(
            view, 'export_gui_to_config'
        ), patch.object(app.isodeceng, 'export_peaks') as export:
            app.batch_process(str(spectrum_path))
        output = Path(directory) / 'sample_unidecfiles'
        assert (output / 'seq.fasta').is_file()
        assert [call.kwargs.get('filename', call.args[1] if len(call.args) > 1 else None)
                for call in export.call_args_list] == [
                    str(output / 'results'), str(output / 'results.tsv')]
        assert os.getcwd() == current_directory
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

    def test_brute_force_uses_prepared_centroids(self):
        from contextlib import redirect_stdout
        from io import StringIO
        from types import SimpleNamespace
        from unittest.mock import Mock, patch

        import numpy as np
        from unidec.IsoDecGUI import IsoDecPres
        from isodec.config import IsoDecConfig

        value = lambda item: SimpleNamespace(GetValue=lambda: item)
        spectrum = np.array([[300.0, 2.0], [300.5, 5.0], [301.0, 3.0]])
        controls = SimpleNamespace(
            ctlavgpeakmasses=value(False), ctlfragmentppm=value("5"),
            ctlfragmentation=value("ETD"),
        )
        view = SimpleNamespace(
            controls=controls, export_gui_to_config=Mock(), SetStatusText=Mock(),
            peakpanel=SimpleNamespace(clear_list=Mock()), clear_all_plots=Mock(),
            plot1=SimpleNamespace(flag=True, data=spectrum, subplot1=SimpleNamespace(
                lines=[object()], texts=[], collections=[], patches=[])),
            plot2=SimpleNamespace(clear_plot=Mock()), clear_fragment_plot=Mock(),
        )
        collection = SimpleNamespace(peaks=[])
        matcher = Mock(return_value=collection)
        pres = IsoDecPres.__new__(IsoDecPres)
        pres.view = view
        pres.eng = SimpleNamespace(
            data=SimpleNamespace(data2=spectrum, rawdata=spectrum + [100, 0]),
        )
        pres.isodeceng = SimpleNamespace(config=IsoDecConfig(), brute_force_pep_match=matcher)
        pres._spectrum_plot_artists = (1, 0, 0, 0)
        pres._sequence_text = Mock(return_value="PEPTIDE")
        pres._save_sequence = Mock()
        pres.translate_config = Mock()
        pres.makeplot1 = Mock()

        output = StringIO()
        with patch("unidec.IsoDecGUI.match_fragments"), redirect_stdout(output):
            pres.on_brute_force_match()

        self.assertIs(matcher.call_args.args[1], spectrum)
        self.assertTrue(matcher.call_args.kwargs["centroided"])
        self.assertIn("Brute Force Match Done. Time:", output.getvalue())
        pres.makeplot1.assert_not_called()
        view.clear_all_plots.assert_not_called()

        mz = np.arange(300.0, 302.0, 0.001)
        intensity = sum(np.exp(-((mz - center) / 0.01) ** 2)
                        for center in (300.3, 300.9, 301.5))
        dense_data = np.column_stack((mz, intensity))
        pres.eng.data.data2 = dense_data
        with patch("unidec.IsoDecGUI.match_fragments"):
            pres.on_brute_force_match()

        self.assertLess(len(matcher.call_args.args[1]), len(dense_data))
        self.assertTrue(matcher.call_args.kwargs["centroided"])
        pres.makeplot1.assert_called_once()

        pres.eng.data.data2 = spectrum
        view.plot1.subplot1.lines.append(object())
        with patch("unidec.IsoDecGUI.match_fragments"):
            pres.on_brute_force_match()
        self.assertEqual(pres.makeplot1.call_count, 2)


if __name__ == "__main__":
    unittest.main()
