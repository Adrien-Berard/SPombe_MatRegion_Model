import random
import numpy as np
import pandas as pd

# Parameters for simulation
chromatine_size = 60
polymerase_count = 0
simulation_steps = 200
adding_position = 15
end_of_replication_position = chromatine_size - 15

vicinity_size = 5
recruitment_probability = 1

change_probability = 0.3
F = change_probability / (1 - change_probability)

# Number of simulations
num_simulations = 1000

# Linear function parameters for adding_poly_proba
slope = 1e-5
intercept = 1e-1

# Polymerase movement probabilities
left_movement_probability = 1/2
right_movement_probability = 1/2

# Set seed for reproducibility
np.random.seed(42)

class Chromatine:
    def __init__(self, histones_count):
        # Initialize chromatine with histones
        self.histones = np.full(histones_count, 'M', dtype='U1')

        # Randomly set approximately 30% of histones to 'U' (unmodified)
        num_unmodified = int(0.3 * histones_count)
        unmodified_positions = np.random.choice(histones_count, size=num_unmodified, replace=False)
        self.histones[unmodified_positions] = 'U'

    def regenerate_histones(self, positions_to_delete):
        unmodified_positions = np.where(self.histones == 'U')[0]

        # Ensure that the size of the sample is not greater than the size of the population
        if len(unmodified_positions) >= len(positions_to_delete):
            selected_positions = np.random.choice(unmodified_positions, size=len(positions_to_delete), replace=False)
            self.histones[selected_positions] = 'A'

    def add_polymerases(self, count, existing_polymerase_positions, adding_position):
        # Add a specified number of new polymerases at non-overlapping positions
        for _ in range(count):
            new_position = adding_position
            while new_position in existing_polymerase_positions:
                new_position += 1
            existing_polymerase_positions.append(new_position)

    def adding_poly_proba(self, adding_position):
        # Linear function of the local density of histones
        # Let's calculate local density as the count of active histones in the vicinity of the polymerase
        start_index = max(0, adding_position - vicinity_size)
        end_index = min(len(self.histones), adding_position + vicinity_size + 1)

        local_density = np.sum(np.isin(self.histones[start_index:end_index], ['M', 'A']))

        # Calculate the probability of adding a new polymerase
        probability = slope * local_density + intercept

        return probability

    def change_next_histones(self, position, p_recruitment, p_change, nth_neighbor=1, vicinity_size=5):
        # Simulate the influence of neighbors on the next histones
        if 1 <= position < len(self.histones) - 1:
            current_histone = self.histones[position]

            # Calculate the influence of vicinity on the recruitment probability
            adjusted_p_recruitment = p_recruitment / (nth_neighbor)  
            if adjusted_p_recruitment > 1:
                adjusted_p_recruitment = 1

            # Probabilistically recruit an enzyme
            if np.random.random() < adjusted_p_recruitment:
                # Apply changes with probability p_change
                if np.random.random() < p_change:
                    nth_position = position + nth_neighbor
                    if nth_position < len(self.histones):
                        nth_histone = self.histones[nth_position]

                        #if current_histone == 'U' and nth_histone == 'A':
                        #    self.histones[nth_position] = 'U'
                        #elif current_histone == 'U' and nth_histone == 'M':
                        #    self.histones[nth_position] = 'U'
                        if current_histone == 'A' and nth_histone == 'U':
                            self.histones[nth_position] = 'A'
                        elif current_histone == 'A' and nth_histone == 'M':
                            self.histones[nth_position] = 'U'
                        elif current_histone == 'M' and nth_histone == 'A':
                            self.histones[nth_position] = 'U'
                        elif current_histone == 'M' and nth_histone == 'U':
                            self.histones[nth_position] = 'M'

class Polymerase:
    def __init__(self, chromatine, position=10, temperature=1.0):
        # Initialize polymerase with a reference to the chromatine, a starting position, and a temperature for movement
        self.chromatine = chromatine
        self.position = position
        self.temperature = temperature

    def delete(self):
        polymerases.remove(self)

    def move(self, chromatine):
        # Define two possible states for movement (left and right)
        states = [0, 1]

        next_position = self.position + 1  # Check only one position ahead

        # Check if there is another polymerase in the next position
        if next_position in existing_polymerase_positions:
            # Do not move if there is another polymerase 1 place after
            return

        probabilities = [left_movement_probability, right_movement_probability]

        # Normalize probabilities to sum to 1
        total_prob = np.sum(probabilities)
        normalized_probabilities = probabilities / total_prob

        # Choose the next position based on normalized probabilities
        self.position = np.random.choice([self.position, next_position], p=normalized_probabilities)

        # Bounding conditions
        if self.position == end_of_replication_position:
            self.delete()

    def change_histones(self, chromatine):
        # Simulate the histone change process by polymerase
        if 0 <= self.position < len(chromatine.histones) and chromatine.histones[self.position] == 'M':
            chromatine.histones[self.position] = 'A'

# Initialize lists to store results of each simulation
results_list = []

# Run simulations
for simulation in range(num_simulations):
    print(simulation)

    # Initialize chromatine and polymerases with a specified temperature
    chromatine = Chromatine(chromatine_size)
    polymerases = [Polymerase(chromatine, temperature=1.0) for _ in range(polymerase_count)]

    # Track existing polymerase positions using a list to avoid duplicates
    existing_polymerase_positions = [polymerase.position for polymerase in polymerases]

    # Lists to store the number of transitions and active histones over time
    transitions_dict = {}

    # Run the simulation
    for frame in range(simulation_steps):
        deleted_positions = []  
        polymerase_positions = []  

        for polymerase in polymerases:
            polymerase.move(chromatine)
            polymerase.change_histones(chromatine)
            polymerase_positions.append(polymerase.position)
            deleted_positions.append(polymerase.position)

        for position in range(1, chromatine_size):
            chromatine.change_next_histones(position, p_recruitment=recruitment_probability, p_change=change_probability, nth_neighbor=np.random.randint(1, chromatine_size), vicinity_size=vicinity_size)

        if np.random.random() < chromatine.adding_poly_proba(adding_position):
            chromatine.add_polymerases(1, existing_polymerase_positions, adding_position)
            new_polymerase_positions = existing_polymerase_positions[-1:]
            new_polymerases = [Polymerase(chromatine, position=pos, temperature=1.0) for pos in new_polymerase_positions]
            polymerases.extend(new_polymerases)

        active_histone_count = np.sum(np.isin(chromatine.histones, ['M', 'A']))

    # Append results to the list
    results_list.append({
        'Simulation': simulation + 1,
        'Transitions': transitions_dict
    })

# Save results to CSV
df_results = pd.DataFrame(results_list)
df_results.to_csv(f'simulation_results_.csv', index=False)
