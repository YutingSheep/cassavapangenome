import pandas as pd
import matplotlib.pyplot as plt
import argparse
def main(file_paths, output_csv_path, output_plot_path):
    file_names = [
        "Haploid",
        "Duplicated",
        "Error",
        "Collapsed"
    ] Load all files into separate dataframes with appropriate column names
    dfs = []
    for file_path, file_name in zip(file_paths, file_names):
        df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["Chromosome", file_name])
        dfs.append(df)
    # Merge all dataframes on the 'Chromosome' column
    plt.figure(figsize=(14, 8))
    ax = combined_df.set_index('Chromosome').plot(kind='bar', colormap='viridis', edgecolor='black', figsize=(14, 8))
    plt.title('Chromosome Statistics', fontsize=16)
    plt.xlabel('Chromosome', fontsize=14)
    plt.ylabel('Values (MB)', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title='Statistics', title_fontsize='13', fontsize='12', loc='upper right')
    # Add MB unit to y-axis labels
    labels = [item.get_text() for item in ax.get_yticklabels()]
    #ax.set_yticklabels([f'{label} MB' for label in labels])

    plt.tight_layout()

    # Save the plot to a file

    plt.tight_layout()

    # Save the plot to a file

    plt.tight_layout()

    # Save the plot to a file

    plt.tight_layout()

    # Save the plot to a file
    plt.savefig(output_plot_path, dpi=300)

    # Show the plot
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process chromosome statistics data.')
    parser.add_argument('files', metavar='F', type=str, nargs=4, help='Paths to the input files')
    parser.add_argument('--output_csv', type=str, default='combined_stat_corrected.csv',
 help='Path to save the combined CSV file')
    parser.add_argument('--output_plot', type=str, default='chromosome_statistics_bar_chart.png', help='Path to save the plot image')

    args = parser.parse_args()
    main(args.files, args.output_csv, args.output_plot)