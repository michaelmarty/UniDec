import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    # Function to read CSV files and plot chromatograms with normalization and baseline correction
    def plot_chromatograms(file_paths, legend_names=None,
                           x_label='Retention Volume (mL)', x_label_fontsize=12, x_label_bold=False,
                           y_label_fontsize=12, y_label_bold=False,
                           legend_fontsize=12, legend_bold=False,
                           line_width=1):
        plt.figure(figsize=(10, 6))

        for idx, file_path in enumerate(file_paths):
            # Read the CSV file, specifying the row containing column names (0-indexed)
            data = pd.read_csv(file_path, header=1)

            # Update these variable names based on the actual column names in your CSV
            volume_column_name = data.columns[0]  # The first column name
            uv_absorbance_column_name = data.columns[1]  # The second column name

            # Assuming the CSV file has the correct columns
            volume = data[volume_column_name]
            uv_absorbance = data[uv_absorbance_column_name]

            # Shift the x-values (volume) of the second dataset (index 1) to the right by 0.30759
            if idx == 1:  # The second file (index 1 in 0-based index)
                shift_value = 0.30759
                volume = volume + shift_value

            # Filter the data between 5 and 20 minutes
            valid_data = (volume >= 5) & (volume <= 20)
            volume_filtered = volume[valid_data]
            uv_absorbance_filtered = uv_absorbance[valid_data]

            # Baseline correction: subtract the minimum value from the absorbance
            baseline_corrected = uv_absorbance_filtered - uv_absorbance_filtered.min()

            # Normalize the baseline-corrected UV absorbance values
            uv_absorbance_normalized = baseline_corrected / baseline_corrected.max()

            # Label each plot with a number or custom name from legend_names
            label = legend_names[idx] if legend_names else f'Chromatogram {idx + 1}'

            # Plot chromatogram with a different colored line for each file
            plt.plot(volume_filtered, uv_absorbance_normalized, label=label, linewidth=line_width)

        # Set the y-axis limit from 0 to 1 for normalized data
        plt.ylim(0, 1)

        # Set the x-axis limit from 5 to 20 minutes
        plt.xlim(5, 20)

        # Configure y-axis ticks to display 0 and 100 as integers
        plt.yticks([0, 1], ['0', '100'], fontsize=y_label_fontsize)

        # Add a bold custom y-axis label centered with the % symbol
        y_label_weight = 'bold' if y_label_bold else 'normal'
        plt.ylabel('%', fontsize=y_label_fontsize, fontweight='bold', labelpad=20)

        # Add x-axis label with adjustable font size and bolding
        x_label_weight = 'bold' if x_label_bold else 'normal'
        plt.xlabel(x_label, fontsize=x_label_fontsize, fontweight=x_label_weight)

        # Add legend with customizable font size and bolding
        legend_weight = 'bold' if legend_bold else 'normal'
        plt.legend(loc='best', prop={'size': legend_fontsize, 'weight': legend_weight})

        # Remove the top and right borders of the plot
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Show the plot
        plt.show()

    # Replace the list with the actual paths to your CSV files
    file_paths = [
        'mewEC_2D_07_2_2024.csv',
        '20231215_HSJ_AqpZWT_E3T0-_EcoliPolarLipid_ND_SEC_Day2.csv',
        'mewEC_MEMB_EX_4A_06_21_2024.csv',
        # Add your actual file paths here
    ]

    # Create full file paths
    base_dir = "D:\\FPLC"
    full_paths = [os.path.join(base_dir, file) for file in file_paths]

    # Provide custom legend names corresponding to each file path
    legend_names = [
        'empty control',
        'AqpZ purified ND',
        'AqpZ membrane ND',
    ]

    # Call the function with customizations
    plot_chromatograms(
        full_paths,
        legend_names,
        x_label='Retention Volume (mL)',
        x_label_fontsize=14, x_label_bold=True,  # Bold x-axis label
        y_label_fontsize=14, y_label_bold=True,  # Bold y-axis label
        legend_fontsize=12, legend_bold=True,   # Bold legend text
        line_width=2                            # Custom line width
    )

main()
