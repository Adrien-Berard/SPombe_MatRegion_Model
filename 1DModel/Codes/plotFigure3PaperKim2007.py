import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming you have a CSV file with the data, replace 'your_csv_file.csv' with your actual CSV file name.
csv_filename = 'fig3_C_coop_counting_lists_dataframe_polymerasecount_0_alpha_0.9871794871794872_F_77_addingpolyprobaintercept_0_addingpolyprobaslope_0.csv'
df = pd.read_csv(csv_filename)

sns.kdeplot(df['Methylated Histone Count'] - df['Acetylated Histone Count'],multiple = 'stack')

plt.xlabel('M-A')
plt.ylabel('P(M-A)')
plt.title("F = 77 coop C")
plt.semilogy()
plt.show()

