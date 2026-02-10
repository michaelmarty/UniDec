import matplotlib.pyplot as plt
import numpy as np
import os


def create_plate_map():
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    cols = [str(i) for i in range(1, 13)]

    # Initialize all wells as empty/white
    well_data = {(r, int(c)): ("", "white") for r in rows for c in cols}
    color_map = {'r': 'red', 'b': 'skyblue', 'g': 'lightgreen', 'y': 'yellow', 'w': 'white'}

    print("--- 96-Well Plate Mapper ---")
    print("Commands: 's'=skip, 'f'=finish/save, 'g'=goto specific well")

    # Iterate through wells but allow jumping
    row_idx, col_idx = 0, 0
    while row_idx < 8:
        r = rows[row_idx]
        c = int(cols[col_idx])
        well_id = f"{r}{c}"

        user_input = input(f"Content for {well_id} (s/f/g/text): ").strip()

        if user_input.lower() == 'f':
            break
        elif user_input.lower() == 'g':
            target = input("Enter well to jump to (e.g., B5): ").strip().upper()
            if len(target) >= 2 and target[0] in rows and target[1:].isdigit():
                row_idx = rows.index(target[0])
                col_idx = int(target[1:]) - 1
                continue
            else:
                print("Invalid well format.")
                continue
        elif user_input.lower() == 's':
            pass  # Just move to next
        else:
            color_choice = input(f"Color for {well_id} (r/g/b/y/w): ").lower()
            well_data[(r, c)] = (user_input, color_map.get(color_choice, 'white'))

        # Increment logic
        col_idx += 1
        if col_idx >= 12:
            col_idx = 0
            row_idx += 1

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.5, 8.5)

    for (r_label, c_num), (text, color) in well_data.items():
        y_pos = 8 - rows.index(r_label)  # Flip Y so A is at the top
        x_pos = c_num

        circle = plt.Circle((x_pos, y_pos), 0.4, color=color, ec='black', lw=1)
        ax.add_patch(circle)
        ax.text(x_pos, y_pos, text, ha='center', va='center', fontsize=7, fontweight='bold')

    # Formatting
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(cols)
    ax.set_yticks(range(1, 9))
    ax.set_yticklabels(reversed(rows))
    ax.set_title(f"Plate Map", fontsize=16)
    plt.gca().set_aspect('equal')

    # Save logic
    outdir = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\96 Well Plate Maps"
    if not os.path.exists(outdir):
        print(f"Warning: Path {outdir} not found. Saving to current directory.")
        outdir = "."

    file_name = input("Enter filename (include .pdf or .png): ")
    plt.savefig(os.path.join(outdir, file_name), dpi=300, bbox_inches='tight')
    print(f"Saved to: {os.path.join(outdir, file_name)}")
    plt.show()


if __name__ == "__main__":
    create_plate_map()