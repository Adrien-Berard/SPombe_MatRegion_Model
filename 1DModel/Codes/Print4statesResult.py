import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# Assuming you have a CSV file with the data, replace 'your_csv_file.csv' with your actual CSV file name.
csv_filename = 'ModellingChromatine/csv/3statesModel_F_100.0.csv'
df = pd.read_csv(csv_filename)
length_chro = 200

# Select the frame you want to visualize
frame_to_show = len(df) - 1

# Subplot for Polymerase Count
plt.subplot(3, 1, 1)
plt.plot(df['Time Steps'][:frame_to_show + 1], 
         df['Polymerase Count'][:frame_to_show + 1], 
         label='Polymerase Count')
plt.title('Polymerase Count')
plt.xlabel('Time Steps in s')
plt.ylabel('Count')
plt.legend()

# Subplot for Histone Counts
plt.subplot(3, 1, 2)
plt.plot(df['Time Steps'][:frame_to_show + 1], 
         df['A Count'][:frame_to_show + 1], 
         label='A', color="b")
plt.plot(df['Time Steps'][:frame_to_show + 1], 
         df['B Count'][:frame_to_show + 1], 
         label='B', color="r")
plt.plot(df['Time Steps'][:frame_to_show + 1], 
         df['C Count'][:frame_to_show + 1], 
         label='C', color="y")
plt.plot(df['Time Steps'][:frame_to_show + 1], 
         df['D Count'][:frame_to_show + 1], 
         label='D', color="g")
plt.title('Histone Counts')
plt.xlabel('Time Steps in s')
plt.ylabel('Count')
plt.legend()


# Adjust layout to prevent subplot overlap
plt.tight_layout()

# Display the plot
plt.show()