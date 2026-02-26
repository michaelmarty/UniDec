# run_plot.py
from cdms_mass_plot_utils import plot_mass_offset_stack

files = [
    r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\sfCherry - done\Long - done\20251208_EC_D1cherryLong_1to90_DMPC_CDMS_2025-12-08-11-04-12_unidecfiles\20251208_EC_D1cherryLong_1to90_DMPC_CDMS_2025-12-08-11-04-12_mass.txt",
    r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\mCherry - done\20251202_EC_D1mCherryREDOcdms_2025-12-02-11-15-03_unidecfiles\20251202_EC_D1mCherryREDOcdms_2025-12-02-11-15-03_mass.txt",
    r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\mBlueberry - done\20251107_EC_D1blueberry_1to80_DMPC_CDMS_2025-11-07-12-07-14_unidecfiles\20251107_EC_D1blueberry_1to80_DMPC_CDMS_2025-11-07-12-07-14_mass.txt",
    r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\Halo-Tag - done\20251204_EC_D1Halo_1to80_DMPC_rep3cdms_2025-12-04-11-13-01_unidecfiles\20251204_EC_D1Halo_1to80_DMPC_rep3cdms_2025-12-04-11-13-01_mass.txt",
    r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\sfGFP - done\E3T0\20251209_EC_E3sfGFP_1to160_DMPC_run_twoCDMS_2025-12-09-01-33-49_unidecfiles\20251209_EC_E3sfGFP_1to160_DMPC_run_twoCDMS_2025-12-09-01-33-49_mass.txt",
]

labels = [
    "MSP1D1-sfCherryL",
    "MSP1D1-mCherry",
    "MSP1D1-mBlueberry",
    "MSP1D1-HaloTag",
    "MSP1E3D1-sfGFP",
]

second_labels = [
    "185 ± 16 DMPC per Nanodisc",
    "154 ± 9 DMPC per Nanodisc",
    "156 ± 7 DMPC per Nanodisc",
    "146 ± 14 DMPC per Nanodisc",
    "360 ± 20 DMPC per Nanodisc",
]

fig, axes = plot_mass_offset_stack(
    datasets=files,
    labels=labels,
    second_labels=second_labels,
    x_range=(0, 500),
    colors=["#CC79A7", "#CC79A7", "#2C7FB8", "#7A7A7A", "#00A9A5"],
    y_offset_step=120,
    y_labels=("0", "%", "100"),
    font_name="Arial",
    font_size=10,
    figsize=(9, 5),
    dpi=300,
    save=None,   # e.g., "output.png"
    show=True
)