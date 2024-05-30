import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the result data from the CSV file
csv_filename = 'averageNEW0.3_m_a_ratio_results_simulationSteps_50000_startingFvalue_0.1_endingFvalue_100.csv'
result_summary_df = pd.read_csv(csv_filename)

# Plot the results
plt.figure(figsize=(10, 6))
plt.semilogx()
plt.plot(result_summary_df['F'], result_summary_df['Average_M_A_ratio'], marker='o', linestyle='-', color='b')
plt.title('Average |M - A| / (M + A) vs. F')
plt.xlabel('F')
plt.ylabel('Average |M - A| / (M + A)')
plt.grid(True)
plt.show()
