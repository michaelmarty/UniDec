from unidec.modules.gui_elements import CDWindow, ChromCDMenu, ChromCD_controls


class CDMainwindow(CDWindow.CDMainwindow):
    """UniChromCD window with chromatography-specific plots, panels, and controls."""

    chrom_mode = True
    controls_class = ChromCD_controls.main_controls
    menu_class = ChromCDMenu.CDMenu
