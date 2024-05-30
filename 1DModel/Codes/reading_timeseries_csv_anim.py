import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Read the CSV file into a DataFrame
csv_filename = 'TimeseriesBurstPERIODChromatin_burstFrequency0.95_burstSize6_ProbaRight0.001_burstDuration2000_InactiveDuration5000.csv'
df = pd.read_csv(csv_filename)

# Parameters for simulation
chromatine_size = 60
polymerase_count = 0

# Function to initialize chromatine
def initialize_chromatine(histones_count):
    histones = np.full(histones_count, 'M', dtype='U1')
    num_unmodified = int(0.5 * histones_count)
    unmodified_positions = np.random.choice(histones_count, size=num_unmodified, replace=False)
    histones[unmodified_positions] = 'U'
    return histones

# Function to visualize chromatine structure
def visualize_chromatine(histones, polymerase_positions=None):
    plt.bar(range(len(histones)), [1] * len(histones),
            color=[('olive' if hist == 'U' else 'red' if hist == 'M' else 'blue' if hist == 'A' else 'green') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")

    if polymerase_positions:
        for position in polymerase_positions:
            plt.arrow(position, -0.2, 0, -0.1, color='black', head_width=0.5, head_length=0.1)

# Function to update the plot in each animation frame
def update(frame):
    plt.clf()  # Clear the previous frame

    # Subplot for Polymerase Count
    plt.subplot(3, 1, 1)
    plt.plot(df['Time Steps'][:frame+1], df['Polymerase Count'][:frame+1], marker='o', label='Polymerase Count')
    plt.title('Polymerase Count')
    plt.xlabel('Time Steps')
    plt.ylabel('Count')
    plt.legend()

    # Subplot for Histone Counts
    plt.subplot(3, 1, 2)
    plt.plot(df['Time Steps'][:frame+1], df['Active Histone Count'][:frame+1], marker='o', label='Active Histone Count')
    plt.plot(df['Time Steps'][:frame+1], df['Acetylated Histone Count'][:frame+1], marker='o', label='Acetylated Histone Count')
    plt.plot(df['Time Steps'][:frame+1], df['Methylated Histone Count'][:frame+1], marker='o', label='Methylated Histone Count')
    plt.plot(df['Time Steps'][:frame+1], df['Unmodified Histone Count'][:frame+1], marker='o', label='Unmodified Histone Count')
    plt.title('Histone Counts')
    plt.xlabel('Time Steps')
    plt.ylabel('Count')
    plt.legend()

    # Subplot for Chromatine Visualization
    plt.subplot(3, 1, 3)
    histones = initialize_chromatine(chromatine_size)
    polymerase_positions = df['Polymerase Count'][frame]
    visualize_chromatine(histones, polymerase_positions)
    plt.title('Chromatine Structure')

# Create an animated plot
fig = plt.figure(figsize=(10, 15))
ani = FuncAnimation(fig, update, frames=len(df), interval=200, repeat=False)

# Show the animation
plt.tight_layout()
plt.show()
