import numpy as np
import pandas as pd

# Parameters for simulation
chromatine_size = 60
polymerase_count = 0
simulation_steps = 50000
adding_position = 25
end_of_replication_position = chromatine_size

# Simulation-specific parameters
start_F = 0.1
end_F = 100
F_values = np.geomspace(start_F,end_F)

# Linear function parameters
slope = 0
intercept = 0

# Polymerase movement probabilities
left_movement_probability = 1/2
right_movement_probability = 1/2

# Set seed for reproducibility
np.random.seed(42)

class Chromatine:
    def __init__(self, histones_count):
        self.histones = np.full(histones_count, 'A', dtype='U1')

        num_unmodified = int(histone_modification_percentage * histones_count)
        unmodified_positions = np.random.choice(histones_count, size=num_unmodified, replace=False)
        self.histones[unmodified_positions] = 'U'


    def noisy_transition(self, position, noisy_transition_probability, noisy_changes):
        current_histone = self.histones[position]
        if np.random.random() < 1/3:
            if current_histone in ['A', 'M']:
                self.histones[position] = 'U'
                noisy_changes += 1
            elif current_histone == 'U':
                self.histones[position] = np.random.choice(['A', 'M'])
                noisy_changes += 1

        return noisy_changes

  
    def change_next_histones(self, position, p_recruitment, p_change, enzyme_changes, nth_neighbor):
        if 1 <= position < len(self.histones) - 1:
            current_histone = self.histones[position]

            nth_position = nth_neighbor
            if nth_position <= len(self.histones):
                nth_histone = self.histones[nth_position]

                if current_histone == 'A' and nth_histone == 'U':
                    self.histones[nth_position] = 'A'
                    enzyme_changes += 1 
                elif current_histone == 'A' and nth_histone == 'M':
                    self.histones[nth_position] = 'U'
                    enzyme_changes += 1 
                elif current_histone == 'M' and nth_histone == 'U':
                    self.histones[nth_position] = 'M'
                    enzyme_changes += 1 
                elif current_histone == 'M' and nth_histone == 'A':
                    self.histones[nth_position] = 'U'
                    enzyme_changes += 1 
        return enzyme_changes

# Create an empty dataframe to store the results for different F values
columns = ['F', 'Average_M_A_ratio']
result_summary_df = pd.DataFrame(columns=columns)

# Set seed for reproducibility
np.random.seed(42)

# Number of simulations for each F value
num_simulations = 3

for F in F_values:
    print(F)
    alpha = F / (1 + F)

    histone_modification_percentage = 0.5
    recruitment_probability = 1
    change_probability = alpha
    regeneration_probability = 0.3
    adding_polymerase_probability = 0.3
    noisy_transition_probability = 1 - alpha
    vicinity_size = 5

    chromatine = Chromatine(chromatine_size)

    m_a_ratios = []

    for simulation in range(num_simulations):
        for frame in range(simulation_steps):
            deleted_positions = []
            polymerase_positions = []
            noisy_changes_count = 0
            enzyme_changes_count = 0

            position = np.random.randint(1, chromatine_size)
            if np.random.random() < alpha:
                enzyme_changes_count = chromatine.change_next_histones(position, p_recruitment=recruitment_probability,
                                                                p_change=change_probability, enzyme_changes=enzyme_changes_count,
                                                                nth_neighbor=np.random.randint(1, chromatine_size))
            else:
                noisy_changes_count = chromatine.noisy_transition(position, noisy_transition_probability, noisy_changes_count)

            active_histone_count = np.sum(np.isin(chromatine.histones, ['M', 'A']))
            acetylated_histone_count = np.sum(chromatine.histones == 'A')
            methylated_histone_count = np.sum(chromatine.histones == 'M')

            m_a_ratio = np.abs(methylated_histone_count - acetylated_histone_count) / active_histone_count
            m_a_ratios.append(m_a_ratio)

        average_m_a_ratio = np.mean(m_a_ratios)
        result_summary_df = pd.concat([result_summary_df, pd.DataFrame([{'F': F, 'Average_M_A_ratio': average_m_a_ratio}])], ignore_index=True)

# Save the dataframe to a CSV file
csv_filename = f'averageNEW0.3_m_a_ratio_results_simulationSteps_{simulation_steps}_startingFvalue_{start_F}_endingFvalue_{end_F}.csv'
result_summary_df.to_csv(csv_filename, index=False)

print("Done")


# 