from chromatogram_plot_utils import plot_sec_stack

files = [
    r"D:\FLND AKTA files\XLSX chromatograms\D1_1to80_DMPC.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to60_DMPC.xlsx",
    r"D:\FLND AKTA files\XLSX chromatograms\D1sfCherrySHORT_1to80_DMPC.xlsx",
    r"D:\FLND AKTA files\New folder\MSP1d1_CHERRYLONG.xlsx",
    r"D:\FLND AKTA files\New folder\mcherry.xlsx",
    r"D:\FLND AKTA files\XLSX chromatograms\d1bLUEBERRY_1TO80_dmpc.xlsx",
    r"D:\FLND AKTA files\New folder\halotag.xlsx",
    r"D:\FLND AKTA files\XLSX chromatograms\E3sfGFP_1to150_DMPC.xlsx",
]

sample_labels = [
    "MSP1D1",
    "MSP1D1-sfGFP",
    "MSP1D1-sfCherryshort",
    "MSP1D1-sfCherryLong",
    "MSP1D1-mCherry",
    "MSP1D1-mBlueberry2",
    "MSP1D1-HaloTag",
    "MSP1E3D1-sfGFP"
]

colors = [
    "#333333", "#00A9A5", "#C76BA0", "#C76BA0",
    "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5",
]

standard_volumes = [13.81, 16.82, 18.09, 19.05, 20.79]

fig, axes = plot_sec_stack(
    files=files,
    sample_labels=sample_labels,
    colors=colors,
    standard_volumes=standard_volumes,
    xlim=(10, 20),
    ylim=(0, 1.0),
    font_name="Arial",
    font_size=10,
    figsize=(8, 10),
    dpi=300,
    save=None,  # e.g., "sec_stack.png"
    show=True
)