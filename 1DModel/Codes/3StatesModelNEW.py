import numpy as np
import pandas as pd
import os

# Parameters for simulation
chromatine_size = 200
polymerase_count = 0
simulation_steps = 1000000
adding_position = 25
end_of_replication_position = chromatine_size - 25

# Simulation-specific parameters
F = 100.0
alpha = F/(1+F)

histone_modification_percentage = 0.5
recruitment_probability = 1
# alpha = 9/10
change_probability = alpha
regeneration_probability = 0.3
adding_polymerase_probability = 0.3
noisy_transition_probability = 1 - alpha
vicinity_size = 5


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
        # Initialize chromatine with histones
        self.histones = np.random.choice(['A', 'B', 'C'], histones_count)

    def noisy_transition(self, position, noisy_transition_probability, noisy_changes):
        if np.random.random() < 1/3:
            if self.histones[position] == 'A':
                self.histones[position] = 'B'
                noisy_changes += 1
            elif self.histones[position] == 'C':
                self.histones[position] = 'B'
                noisy_changes += 1
            elif self.histones[position] == 'B':
                if np.random.random() < 1/2:
                    self.histones[position] = 'A'
                    noisy_changes += 1
                else:
                    self.histones[position] = 'C'
                    noisy_changes += 1

        return noisy_changes

    # def add_polymerases(self, count, existing_polymerase_positions, adding_position):
    #     for _ in range(count):
    #         new_position = adding_position
    #         while new_position in existing_polymerase_positions:
    #             new_position += 1
    #         if new_position < end_of_replication_position:
    #             existing_polymerase_positions.append(new_position)

    # def adding_poly_proba(self, adding_position):
    #     start_index = max(0, adding_position - vicinity_size)
    #     end_index = min(len(self.histones), adding_position + vicinity_size + 1)
    #     local_density = np.sum(np.isin(self.histones[start_index:end_index], ['M', 'A']))
    #     probability = slope * local_density + intercept
    #     return probability

    def change_next_histones(self, position, p_recruitment, p_change, enzyme_changes, nth_neighbor):
        if 1 <= position < len(self.histones) - 1:
            current_histone = self.histones[position]
            # adjusted_p_recruitment = p_recruitment
            # if adjusted_p_recruitment > 1:
            #     adjusted_p_recruitment = 1

            # if np.random.random() < adjusted_p_recruitment:
                # if np.random.random() < p_change:
            nth_position = position + nth_neighbor
            if nth_position < len(self.histones):
                nth_histone = self.histones[nth_position]

                if current_histone == 'A' and nth_histone == 'B':
                    self.histones[nth_position] = 'A'
                    enzyme_changes += 1 
                elif current_histone == 'A' and nth_histone == 'C':
                    self.histones[nth_position] = 'B'
                    enzyme_changes += 1 
                    enzyme_changes += 1 
                elif current_histone == 'B' and nth_histone == 'A':
                    self.histones[nth_position] = 'B'
                    enzyme_changes += 1 
                elif current_histone == 'B' and nth_histone == 'C':
                    self.histones[nth_position] = 'B'
                    enzyme_changes += 1 
                    enzyme_changes += 1 
                elif current_histone == 'C' and nth_histone == 'A':
                    self.histones[nth_position] = 'B'
                    enzyme_changes += 1
                elif current_histone == 'C' and nth_histone == 'B':
                    self.histones[nth_position] = 'C'
                    enzyme_changes += 1

        return enzyme_changes

# class Polymerase:
#     def __init__(self, chromatine, position=10, temperature=1.0):
#         self.chromatine = chromatine
#         self.position = position
#         self.temperature = temperature

#     def delete(self):
#         polymerases.remove(self)

#     def move(self, chromatine):
#         states = [0, 1]
#         next_position = self.position + 1

#         if next_position in existing_polymerase_positions:
#             return

#         probabilities = [left_movement_probability, right_movement_probability]
#         total_prob = np.sum(probabilities)
#         normalized_probabilities = probabilities / total_prob

#         self.position = np.random.choice([self.position, next_position], p=normalized_probabilities)

#         if self.position == end_of_replication_position:
#             self.delete()

#     def change_histones(self, chromatine):
#         if 0 <= self.position < len(chromatine.histones) and chromatine.histones[self.position] == 'U':
#             chromatine.histones[self.position] = 'A'
#         elif 0 <= self.position < len(chromatine.histones) and chromatine.histones[self.position] == 'M':
#             chromatine.histones[self.position] = 'U'

# Initialize chromatine and polymerases with a specified temperature
chromatine = Chromatine(chromatine_size)
# polymerases = [Polymerase(chromatine, temperature=1.0) for _ in range(polymerase_count)]

# Track existing polymerase positions using a list to avoid duplicates
# existing_polymerase_positions = [polymerase.position for polymerase in polymerases]

# Create an empty dataframe to store the counting lists
columns = ['Time Steps', 'Polymerase Count', 'A Count', 'B Count','C Count','D Count', 'Noisy Changes Count', 'Enzyme Changes Count']
result_df = pd.DataFrame(columns=columns)

# Simulation loop
for frame in range(simulation_steps):

    deleted_positions = []  # Keep track of deleted positions for regeneration
    polymerase_positions = []  # Clear the polymerase_positions list
    noisy_changes_count = 0
    enzyme_changes_count = 0

    # for polymerase in polymerases:
    #     polymerase.move(chromatine)
    #     polymerase.change_histones(chromatine)
    #     polymerase_positions.append(polymerase.position)  # Append the current position
    #     deleted_positions.append(polymerase.position)

    # Change the next histones based on the influence of first neighbors
    position = np.random.randint(1, chromatine_size)
    # Use p_recruitment and p_change probabilities with decreasing probability with vicinity
    if np.random.random() < alpha:
        enzyme_changes_count = chromatine.change_next_histones(position, p_recruitment=recruitment_probability,
                                                          p_change=change_probability, enzyme_changes=enzyme_changes_count,
                                                          nth_neighbor=np.random.randint(1, chromatine_size))
    else:
        noisy_changes_count = chromatine.noisy_transition(position, noisy_transition_probability, noisy_changes_count)

    # Regenerate histones at unmodified positions
    # if np.random.random() < regeneration_probability:
    #   chromatine.regenerate_histones(deleted_positions)

    # Randomly add new polymerase at the beginning of the chromatine with a certain probability
    # if np.random.random() < chromatine.adding_poly_proba(adding_position):
    #     # Add new polymerases with non-overlapping random positions
    #     chromatine.add_polymerases(1, existing_polymerase_positions, adding_position)
    #     new_polymerase_positions = existing_polymerase_positions[-1:]
    #     new_polymerases = [Polymerase(chromatine, position=pos, temperature=1.0) for pos in new_polymerase_positions]
    #     polymerases.extend(new_polymerases)

    # Update the number of polymerases and active histones lists
    # polymerase_count_over_time = len(polymerases)
    A_count = np.sum(chromatine.histones == 'A')
    B_count = np.sum(chromatine.histones == 'B')
    C_count = np.sum(chromatine.histones == 'C')
    D_count = np.sum(chromatine.histones == 'D')


    if frame%100 == 0:
        print(frame)
        # Append data to the dataframe
        result_df = pd.concat([result_df, pd.DataFrame([{'Time Steps': frame + 1,
                                                    'Polymerase Count': 0,
                                                    'A Count': A_count,
                                                    'B Count': B_count,
                                                    'C Count': C_count,
                                                    'D Count': D_count,
                                                    'Noisy Changes Count': noisy_changes_count,
                                                    'Enzyme Changes Count': enzyme_changes_count}])], ignore_index=True)

print("Done")

# Save the dataframe to a CSV file

current_directory = os.getcwd()
csv_filename = f'3statesModel_F_{F}.csv'
result_df.to_csv(csv_filename, index=False)