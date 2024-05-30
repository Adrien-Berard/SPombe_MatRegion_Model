import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Assuming you have a CSV file with the data, replace 'your_csv_file.csv' with your actual CSV file name.
csv_filename = 'ScanCenH_Proba_0.8_polymerasecount_0_F_10_newpolyproba_0.01.csv'
df = pd.read_csv(csv_filename)

# Define the number of simulation steps and frames to show
simulation_steps = 100000
frames_to_show = simulation_steps // 200

# Extract and reshape 'A in gene' column
As = df['A in gene'].to_numpy()
As_reshaped = As.reshape(frames_to_show, len(As) // frames_to_show)
meanA, stdA = As_reshaped.mean(axis=1), As_reshaped.std(axis=1)

burst = np.linspace(0.1, 1, 40)
cenH = range(15, 31)

# Use np.repeat and np.tile to ensure consistent lengths
cenHsize_col = np.repeat(list(cenH), len(burst))
burst_col = np.tile(burst, len(cenH))
mean_col = np.repeat(meanA, len(cenH) * len(burst) // len(meanA))
stdA_col = np.repeat(stdA, len(cenH) * len(burst) // len(stdA))

# Ensure that all arrays have the same length
length_check = len(cenHsize_col) == len(burst_col) == len(mean_col) == len(stdA_col)
if not length_check:
    raise ValueError("All arrays must be of the same length")

# Create a scatter plot using seaborn
g = sns.relplot(
    x=cenHsize_col, y=burst_col,
    hue=mean_col, size=stdA_col,
    palette='viridis',
    height=6,
)

# Set log scale for both x and y axes
g.set(xscale="log", yscale="log")

# Add colorbar
cbar = plt.colorbar(g.ax.collections[0])
cbar.set_label('Mean A in gene')

# Set axis labels and title
plt.xlabel('Burst Frequency (log scale)')
plt.ylabel('CenHsize (log scale)')
plt.title('Scatter Plot of Burst Frequency vs. CenHsize')

# Show the plot
plt.show()
