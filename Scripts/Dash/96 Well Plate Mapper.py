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


def create_summary_page(pdf, exp_details):
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    title = f"Plate Mapping Session: {exp_details['name']}"

    plt.text(0.5, 0.95, title, ha='center', va='center', fontsize=16, fontweight='bold')
    plt.text(0.5, 0.91, f"Generated on: {timestamp}", ha='center', va='center', fontsize=10)

    headers = ["Exp Type", "Core/Time", "Plate", "Flow Well", "Elute Well"]
    table_data = [headers] + exp_details['log']

    table = ax.table(cellText=table_data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)

    pdf.savefig(fig)
    plt.close(fig)


def plot_plate(plate_num, data, rows):
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.5, 8.5)

    for r_idx, r_label in enumerate(rows):
        for c_idx in range(1, 13):
            well_info = data.get((plate_num, r_label, c_idx), ("", "white"))
            y_pos = 8 - r_idx

            # Using the stored color from the data dictionary
            circle = plt.Circle((c_idx, y_pos), 0.4, color=well_info[1], ec='black', lw=1)
            ax.add_patch(circle)
            if well_info[0]:
                ax.text(c_idx, y_pos, well_info[0], ha='center', va='center',
                        fontsize=5, fontweight='bold', rotation=45)

    ax.set_xticks(range(1, 13))
    ax.set_yticks(range(1, 9))
    ax.set_yticklabels(reversed(rows))
    ax.set_title(f"Plate Map - Plate {plate_num + 1}", fontsize=14, pad=20)
    plt.gca().set_aspect('equal')
    return fig


def create_plate_map():
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    all_plates_data = {}
    exp_details = {'name': '', 'log': []}

    print("--- 96-Well Session Manager ---")
    start_well = input("Enter starting well for the whole session (e.g., A1): ").upper()
    p_idx, r_idx = 0, rows.index(start_well[0])
    c_idx = int(start_well[1:]) - 1

    session_active = True
    while session_active:
        print(f"\n--- New Block ---")
        print(f"Current Position: Plate {p_idx + 1}, Well {rows[r_idx]}{c_idx + 1}")
        mode_choice = input("Select Mode (1: Manual, 2: Time): ")

        # Define color pairs for this block
        if mode_choice == '1':
            colors = ("orangered", "lime")
            mode_label = "Manual"
        else:
            colors = ("violet", "aquamarine")
            mode_label = "Time"

        if mode_choice == '1':
            while True:
                core = input("Enter Core Name (or 'done' to finish this block): ")
                if core.lower() == 'done': break

                # Flow Well
                f_w = f"{rows[r_idx]}{c_idx + 1}"
                all_plates_data[(p_idx, rows[r_idx], c_idx + 1)] = (f"flow_{core}", colors[0])
                p_idx, r_idx, c_idx = get_next_well(p_idx, r_idx, c_idx)

                # Elute Well
                e_w = f"{rows[r_idx]}{c_idx + 1}"
                all_plates_data[(p_idx, rows[r_idx], c_idx + 1)] = (f"elute_{core}", colors[1])

                exp_details['log'].append([mode_label, core, p_idx + 1, f_w, e_w])
                p_idx, r_idx, c_idx = get_next_well(p_idx, r_idx, c_idx)

        elif mode_choice == '2':
            core = input("Enter Core Name: ")
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

        more = input("\nAdd another experiment block to this session? (y/n): ").lower()
        if more != 'y':
            session_active = False

    # Save Logic
    outdir = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\96 Well Plate Maps"
    if not os.path.exists(outdir): os.makedirs(outdir, exist_ok=True)

    fname = input("\nEnter filename (e.g., 2026-02-10_Combined_Exp.pdf): ")
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