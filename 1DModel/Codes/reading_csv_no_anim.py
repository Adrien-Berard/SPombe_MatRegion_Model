import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# Assuming you have a CSV file with the data, replace 'your_csv_file.csv' with your actual CSV file name.
csv_filename = 'Lotka-VolteraAnalysis.csv'
df = pd.read_csv(csv_filename)
length_chro = 198

# Select the frame you want to visualize
frame_to_show = len(df) - 1

# Subplot for Polymerase Count
plt.subplot(3, 1, 1)
plt.plot(df['Real Time'][:frame_to_show + 1], 
         df['Polymerase Count'][:frame_to_show + 1], 
         label='Polymerase Count')
plt.title('Polymerase Count')
plt.xlabel('Real Time in s')
plt.ylabel('Count')
plt.legend()

# Subplot for Histone Counts
plt.subplot(3, 1, 2)
plt.plot(df['Real Time'][:frame_to_show + 1], 
         df['Acetylated Histone Count'][:frame_to_show + 1], 
         label='Acetylated', color="b")
plt.plot(df['Real Time'][:frame_to_show + 1], 
         df['Methylated Histone Count'][:frame_to_show + 1], 
         label='Methylated', color="r")
plt.plot(df['Real Time'][:frame_to_show + 1], 
         df['Unmodified Histone Count'][:frame_to_show + 1], 
         label='Unmodified', color="y")
plt.title('Histone Counts')
plt.xlabel('Real Time in s')
plt.ylabel('Count')
plt.legend()

# Subplot for Chromatin State Colormap
plt.subplot(3, 1, 3)
chromatin_states = df['Chromatine Array'][:frame_to_show + 1]

# Convert chromatin states to a numerical representation for colormap
state_to_number = {'U': 0, 'M': 1, 'A': 2}

# Create an array to store each row as a separate column
chromatin_numeric = np.zeros((length_chro, len(chromatin_states)))

for col, chromatin_state in enumerate(chromatin_states):
    # Extract only letters
    letters_only = ''.join(filter(str.isalpha, chromatin_state))
    # Convert chromatin state to a numerical representation for colormap
    chromatin_numeric[:, col] = [state_to_number[char] for char in letters_only]

# Create a colormap with distinct colors for each state
cmap = ListedColormap(['y', 'r', 'b'])

plt.imshow(chromatin_numeric, cmap=cmap, aspect='auto', interpolation='none')
plt.title('Chromatin State Colormap')
plt.xlabel('Not real Time')
plt.ylabel('Chromatin Position')
plt.yticks([64, 94, 131, 141])  # Adjust y-axis ticks for better visualization

# Adjust layout to prevent subplot overlap
plt.tight_layout()

# Display the plot
plt.show()