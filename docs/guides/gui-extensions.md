# Extend the UniDec GUI

These are **Tier 2** Python extension points for contributors building a UniDec desktop tool. Use a display-capable Python environment with the UniDec `gui` extra installed. Create windows and handle wx events on the GUI thread. Use the [headless engines or sessions](../api/index.md) for calculations and automation that do not need a window.

## Launcher to engine

| Layer | Entry point and lifecycle |
| --- | --- |
| Launcher | The `gunidec` command and `python -m unidec.Launcher` call `unidec.Launcher.run_launcher`. The launcher selects a built-in presenter from the command line or a launcher action. A custom presenter can use its own Python entry point and the same `start()` lifecycle. |
| Presenter | A conventional extension can subclass `unidec.GUniDec.UniDecApp`. Its base `unidec.modules.unidec_presbase.UniDecPres` creates the `wx.App`. The conventional presenter's `init()` creates the engine first, then the view, wires events, and reads startup state. Pass `ignore_args=True` for programmatic construction. |
| View | `unidec.modules.mainwindow.Mainwindow` is the conventional wx frame. It holds the presenter in `pres` and receives the engine's config object. Controls, menus, plots, and the peak panel live on the frame. Bind commands to the frame and route analysis actions through the presenter. |
| Engine | The presenter owns the engine as `eng`. The engine owns configuration, input and result arrays, and peaks. It runs processing and native deconvolution. Keep reusable calculations here or in a headless module so they can also serve scripts and tests. |

`start()` sets the top window, shows it, and enters the wx event loop. A presenter already creates its own `wx.App`, so construct it directly in the standalone process responsible for that application. Attach tools and dialogs to that frame and use its event loop. Call `wx.CallAfter` to schedule any worker's UI updates on the GUI thread. Keep bindings and resources scoped to the view's lifetime.

## Hooks to use

- **Add a GUI command.** Call `super().init(...)` in the presenter subclass, then install the command on `self.view`. Use a wx event handler on the presenter to read current engine state or call an existing workflow action. The example binds a menu item this way. If the presenter supports changing plot layouts, `remake_mainwindow` replaces its frame. Reattach custom commands after calling the base remake method.
- **Keep the controls and engine in sync.** `presenter.import_config()` pushes engine settings into controls. `presenter.export_config()` pulls edits into the engine. Passing a filename invokes configuration I/O, and HDF5-based presenters can persist to their open project. The standard workflow handlers call these methods around analysis. Reuse `on_open_file`, `on_dataprep_button`, `on_unidec_button`, and `on_pick_peaks` when an action needs their normal status, plot, and peak-panel updates.
- **Update custom results after a stage.** The conventional `on_unidec_button` calls `after_unidec_run()` after a successful engine run has imported its output. The conventional `on_pick_peaks` calls `after_pick_peaks()` for interactive display updates before scoring. Override the corresponding method and call `super()` to keep built-in updates. Code that needs final scores runs after the scoring action or reads them after the peak-picking handler returns.
- **Add plots in a derived view.** The conventional `Mainwindow` calls `setup_mode_tabbed_plots(plotwindow, figsize)` to create plot widgets parented by its notebook or notebook pages. `add_mode_tabs(plotwindow, section)` inserts pages among the shared tabs at the `"mz"`, `"ccs"`, and `"mass"` locations. In the scrolling layout, use `setup_mode_scrolled_plots(plotwindow, figsize)` and `layout_scrolled_plots()`. The layout must place shared and custom plots. After creating the active layout, `extend_mode_plot_metadata()` can append custom plot widgets and their export names to `plots` and `plotnames` in matching order. The built-in `unidec.modules.IM_mainwindow.IMMainwindow` implements all five hooks.

Match the presenter, view, and engine to the same acquisition mode. For example, the IM-MS presenter creates `UniDecIM` and `IMMainwindow` together. Configuration synchronization preserves the mode selected when the application was built. View controls, internal menu identifiers, and PubSub topics can change with the implementation. For an external tool, anchor custom commands to the documented presenter and layout hooks.

## Example: add a spectrum summary command

This standalone example adds a menu to the conventional presenter and verifies that a real wx menu event reads the opened spectrum from the engine. It uses a temporary workspace for the spectrum, default configuration, and recent-file list. Run it in its own Python process with a graphical display. The verification dispatches the event directly and cleans up the frame, so it needs neither human interaction nor a native deconvolution run.

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import wx

from unidec.GUniDec import UniDecApp
from unidec.modules.unidecstructure import UniDecConfig


class SpectrumSummaryApp(UniDecApp):
    def init(self, *args, **kwargs):
        super().init(*args, **kwargs)
        self.install_summary_menu()

    def install_summary_menu(self):
        menu = wx.Menu()
        self.summary_item = menu.Append(wx.ID_ANY, "Summarize spectrum")
        self.view.GetMenuBar().Append(menu, "Spectrum Summary")
        self.view.Bind(wx.EVT_MENU, self.on_summary, self.summary_item)

    def on_summary(self, event):
        count = len(self.eng.data.rawdata)
        self.view.SetStatusText(f"Spectrum points: {count}", number=5)


with TemporaryDirectory() as directory:
    root = Path(directory)
    source = root / "demo.txt"
    np.savetxt(source, [[500, 1], [501, 3], [502, 2]])

    initialize_paths = UniDecConfig.initialize_system_paths

    def initialize_demo_paths(config):
        initialize_paths(config)
        config.defaultconfig = str(root / "default_conf.dat")
        config.recentfile = str(root / "recent.txt")
        # The presenter falls back to the package list if this file is absent.
        Path(config.recentfile).touch()

    with patch.object(UniDecConfig, "initialize_system_paths", initialize_demo_paths), \
            patch.dict("os.environ", {"UNIDEC_GUI_BRIDGE_FILE": ""}):
        app = SpectrumSummaryApp(ignore_args=True)

    try:
        assert Path(app.eng.config.recentfile) == root / "recent.txt"
        app.on_open_file(str(source))
        assert app.view.config is app.eng.config

        item_id = app.summary_item.GetId()
        event = wx.CommandEvent(wx.EVT_MENU.typeId, item_id)
        assert app.view.GetEventHandler().ProcessEvent(event)
        assert app.view.GetStatusBar().GetStatusText(5) == "Spectrum points: 3"
    finally:
        app.view.Destroy()
        app.wx_app.Yield()
```

In an interactive entry point, construct the subclass and call its inherited `start()`. For a test, use an isolated subprocess so wx state and PubSub listeners are released when the process exits. Dialogs and frames created by an extension should be destroyed when their owning view closes.
