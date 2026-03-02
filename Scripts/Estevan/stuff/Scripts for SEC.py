import pandas as pd
import matplotlib.pyplot as plt
import os



def main():

    # Function to read CSV files and plot chromatograms with custom legend names and axis labels
    def plot_chromatograms(file_paths: object, legend_names: object, x_label: object = None, y_label: object = None) -> object:
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
            # Use the provided legend name for each plot
            label = legend_names[idx]
            # Plot chromatogram with a different colored line for each file
            # Bold the 6th file plot line
            if idx == 5:  # Remember that the index is 0-based, so 5 is the 6th file
                plt.plot(volume, uv_absorbance, label=label, linewidth=1)
            else:
                plt.plot(volume, uv_absorbance, label=label, linewidth=1)
        # Set the y-axis limit from 0 to 500
        plt.ylim(0, 300)
        # Add titles and labels
        plt.title('SEC profile of nanodiscs')
        plt.xlabel(x_label if x_label else volume_column_name)
        plt.ylabel(y_label if y_label else uv_absorbance_column_name)
        # Add legend to the right outside of the plot area
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Show the plot
        plt.show()


    # Replace the list with the actual paths to your CSV files
    file_paths = [
        'mewEC_2D_07_2_2024.csv',
        '20231215_HSJ_AqpZWT_E3T0-_EcoliPolarLipid_ND_SEC_Day2.csv',
        'mewEC_MEMB_EX_4A_06_21_2024.csv',
        # Add your file paths here
    ]

    ls = []
    for file in file_paths:

        curr = "D:\FPLC"
        full_path = os.path.join(curr, file)  # Join the path
        ls.append(full_path)  # Append the joined path



    # Provide custom legend names corresponding to each file path
    legend_names = [
        'empty control',
        'AqpZ purified ND',
        'AqpZ membrane ND',

        # Add your custom names here
    ]

    # Specify custom axis labels if desired
    custom_x_label = 'Retention Volume (mL)'
    custom_y_label = 'UV Absorbance (mAU)'

    plot_chromatograms(ls, legend_names, x_label=custom_x_label, y_label=custom_y_label)


main()