import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from datetime import datetime


def get_next_well(p_idx, r_idx, c_idx):
    c_idx += 1
    if c_idx >= 12:
        c_idx = 0
        r_idx += 1
    if r_idx >= 8:
        r_idx = 0
        p_idx += 1
    return p_idx, r_idx, c_idx


def parse_well(well_str, rows):
    """Helper to convert A1 style input to indices."""
    well_str = well_str.upper().strip()
    try:
        r = rows.index(well_str[0])
        c = int(well_str[1:]) - 1
        return r, c
    except (ValueError, IndexError):
        return None


def create_summary_page(pdf, exp_details):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis('off')
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    plt.text(0.5, 0.95, f"Plate Mapping Session: {exp_details['name']}", ha='center', fontsize=16, fontweight='bold')
    plt.text(0.5, 0.91, f"Generated on: {timestamp}", ha='center', fontsize=10)

    headers = ["Exp Type", "Sample Label", "Plate", "Flow Well", "Elute Well"]
    table_data = [headers] + exp_details['log']
    table = ax.table(cellText=table_data, loc='center', cellLoc='center', colWidths=[0.15, 0.35, 0.1, 0.2, 0.2])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)
    pdf.savefig(fig)
    plt.close(fig)


def plot_plate(plate_num, data, rows):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.5, 8.5)
    for r_idx, r_label in enumerate(rows):
        for c_idx in range(1, 13):
            well_info = data.get((plate_num, r_label, c_idx), ("", "white"))
            y_pos = 8 - r_idx
            circle = plt.Circle((c_idx, y_pos), 0.4, color=well_info[1], ec='black', lw=0.8)
            ax.add_patch(circle)
            if well_info[0]:
                ax.text(c_idx, y_pos, well_info[0], ha='center', va='center', fontsize=5, fontweight='bold',
                        rotation=45)
    ax.set_xticks(range(1, 13))
    ax.set_yticks(range(1, 9))
    ax.set_yticklabels(reversed(rows))
    ax.set_title(f"Plate Map - Plate {plate_num + 1}", fontsize=14, pad=20)
    plt.gca().set_aspect('equal')
    return fig


def create_plate_map():
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    all_plates_data = {}
    exp_details = {'log': []}
    p_idx, r_idx, c_idx = 0, 0, 0  # Initialize at start

    print("--- 96-Well Session Manager ---")

    session_active = True
    while session_active:
        current_well_id = f"{rows[r_idx]}{c_idx + 1}"
        print(f"\n--- NEW EXPERIMENT BLOCK ---")
        print(f"Current Position: Plate {p_idx + 1}, Well {current_well_id}")

        jump = input(f"Press Enter to start at {current_well_id} or type new well (e.g. C1): ").strip().upper()
        if jump:
            coords = parse_well(jump, rows)
            if coords:
                r_idx, c_idx = coords
                print(f"Jumped to {jump}")
            else:
                print(f"Invalid well. Staying at {current_well_id}")

        mode_choice = input("Select Mode: 1 (Manual [OrangeRed/Lime]), 2 (Time [Violet/Aqua]): ").strip()
        colors = ("orangered", "lime") if mode_choice == '1' else ("violet", "aquamarine")
        mode_label = "Manual" if mode_choice == '1' else "Time"

        if mode_choice == '1':
            while True:
                core = input("Enter Core Name (or 'done' to finish block): ").strip()
                if core.lower() == 'done': break
                f_w = f"{rows[r_idx]}{c_idx + 1}"
                all_plates_data[(p_idx, rows[r_idx], c_idx + 1)] = (f"flow_{core}", colors[0])
                p_idx, r_idx, c_idx = get_next_well(p_idx, r_idx, c_idx)
                e_w = f"{rows[r_idx]}{c_idx + 1}"
                all_plates_data[(p_idx, rows[r_idx], c_idx + 1)] = (f"elute_{core}", colors[1])
                exp_details['log'].append([mode_label, core, p_idx + 1, f_w, e_w])
                p_idx, r_idx, c_idx = get_next_well(p_idx, r_idx, c_idx)

        elif mode_choice == '2':
            core = input("Enter Core Name: ").strip()
            end_t = int(input("Enter end time (e.g., 12): "))
            for t in range(end_t + 1):
                t_str = f"t{t:02d}"
                f_w = f"{rows[r_idx]}{c_idx + 1}"
                all_plates_data[(p_idx, rows[r_idx], c_idx + 1)] = (f"flow_{core}_{t_str}", colors[0])
                p_idx, r_idx, c_idx = get_next_well(p_idx, r_idx, c_idx)
                e_w = f"{rows[r_idx]}{c_idx + 1}"
                all_plates_data[(p_idx, rows[r_idx], c_idx + 1)] = (f"elute_{core}_{t_str}", colors[1])
                exp_details['log'].append([mode_label, f"{core}_{t_str}", p_idx + 1, f_w, e_w])
                p_idx, r_idx, c_idx = get_next_well(p_idx, r_idx, c_idx)

        more = input("\nAdd another experiment block? (y/n): ").strip().lower()
        if more not in ['y', 'yes']:
            session_active = False

    # Save Logic
    outdir = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\96 Well Plate Maps"
    if not os.path.exists(outdir): os.makedirs(outdir, exist_ok=True)
    fname = input("\nEnter filename (e.g. 2026-02-11_Exp.pdf): ").strip()
    if not fname.lower().endswith('.pdf'): fname += ".pdf"
    exp_details['name'] = fname
    full_path = os.path.join(outdir, fname)

    with PdfPages(full_path) as pdf:
        create_summary_page(pdf, exp_details)
        for i in range(p_idx + 1):
            fig = plot_plate(i, all_plates_data, rows)
            pdf.savefig(fig)
            plt.close(fig)
    print(f"\nSaved {p_idx + 1} plate(s) to: {full_path}")


if __name__ == "__main__":
    create_plate_map()