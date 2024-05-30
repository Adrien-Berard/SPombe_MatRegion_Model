import numpy as np
import pandas as pd
import os

# Parameters for simulation
chromatine_size = 198
polymerase_count_init = 0
polymerase_count = polymerase_count_init 
simulation_steps = 100000
adding_position = 131
end_of_replication_position = adding_position + 7

gene_position = np.arange(adding_position, end_of_replication_position)

# Simulation-specific parameters
F = 10
alpha = F/(1+F)

histone_modification_percentage = 0.5
recruitment_probability = 1

change_probability = alpha
noisy_transition_probability = 1 - alpha
CenHSart = 65

CenHsize = 30

CenH_positions = np.arange(CenHSart,CenHSart + CenHsize) 
MCenHProb = 0.8

new_poly_probability = 1

# Set seed for reproducibility
np.random.seed(42)

class Chromatine:
    def __init__(self, histones_count,CenH_positions):
        # Initialize chromatine with histones
        self.histones = np.full(histones_count, 'A')
        # full A to start

        # Randomly set approximately histone_modification_percentage of histones to 'U' (unmodified)
        num_unmodified = int(histone_modification_percentage * histones_count)
        unmodified_positions = np.random.choice(histones_count, size=num_unmodified, replace=False)
        self.histones[unmodified_positions] = 'U'
    
        # cenH region
        self.histones[CenH_positions] = 'M'

    def noisy_transition(self, position,CenH_positions, noisy_transition_probability, noisy_changes):
        if position not in CenH_positions:
            if self.histones[position] == 'A':
                self.histones[position] = 'U'
                noisy_changes += 1
            elif self.histones[position] == 'M':
                self.histones[position] = 'U'
                noisy_changes += 1
            elif self.histones[position] == 'U':
                if np.random.random() < 1/2:
                    self.histones[position] = 'A'
                    noisy_changes += 1
                else:
                    self.histones[position] = 'M'
                    noisy_changes += 1
        elif self.histones[position] == 'M' and np.random.random() < 1- MCenHProb:
            self.histones[position] = 'U'
            noisy_changes += 1
        elif self.histones[position] == 'U':
            self.histones[position] = 'M'
            noisy_changes += 1

        return noisy_changes

    def add_polymerases(self, count, adding_position,existing_polymerase_positions):
        for _ in range(count):
            if adding_position not in existing_polymerase_positions:
                existing_polymerase_positions.append(adding_position)
        return existing_polymerase_positions

    def BurstPoly(self,adding_position, num_poly_burst, existing_polymerase_positions):
        Burst_adding_position = adding_position
        for poly in range(num_poly_burst):
            existing_polymerase_positions = chromatine.add_polymerases(1,Burst_adding_position,existing_polymerase_positions)
            Burst_adding_position += 1
        return existing_polymerase_positions

    def adding_poly_proba(self, adding_position, existing_polymerase_positions):
        new_position = adding_position
        probability = new_poly_probability 
        # TO CHANGE AFTERWARDS
        if (new_position in existing_polymerase_positions or new_position >= end_of_replication_position or 
            new_position <= end_of_replication_position) and self.histones[new_position-1] == 'M' and self.histones[new_position-2] == 'M' and self.histones[new_position - 3] == 'M':
                   # Can't bind if 'M-M-M' 
                   # PROBLEMS WITH THE NESTED IF OR IF AND IF
                   probability = 0
        return probability

    def change_next_histones(self,position ,CenH_positions, p_recruitment, p_change, enzyme_changes, nth_neighbor):
        if 1 <= position < len(self.histones) - 1:
            current_histone = self.histones[position]
            nth_position =  nth_neighbor
            nth_histone = self.histones[nth_position]
            if nth_position not in CenH_positions:
                if nth_position < len(self.histones):
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
            elif current_histone == 'M' and nth_histone == 'U' :
                self.histones[nth_position] = 'M'
                enzyme_changes += 1 
            elif current_histone == 'U' and nth_histone == 'M' and np.random.random() < 1- MCenHProb: 
                self.histones[nth_position] = 'U'
                enzyme_changes += 1
        return enzyme_changes

    def CenHRegion(self,CenH_positions, cenHStart, McenHDensity):
        num_no_M = int((1 - McenHDensity) * CenHsize)
        
        unmodified_positions = np.random.choice(CenHsize, size=num_no_M, replace=False)
        
        self.histones[CenH_positions] = 'M'
        
        for position in unmodified_positions:
            self.histones[cenHStart + position] = 'U'
            np.delete(CenH_positions,position)
        
        return CenH_positions

    def ChromatineVisualisation(self):
        return self.histones
    
class Polymerase:
    def __init__(self, chromatine, position=adding_position):
        self.chromatine = chromatine
        self.position = position

    def delete(self,position,existing_polymerase_positions):
        polymerases.remove(self)
        existing_polymerase_positions.remove(position)
        return existing_polymerase_positions
    
    def move(self, chromatine, existing_polymerase_positions):
        states = [0, 1]
        previous_position = self.position
        next_position = self.position + 1

        if next_position not in existing_polymerase_positions:
            probabilities = [left_movement_probability, right_movement_probability]
            total_prob = np.sum(probabilities)
            normalized_probabilities = probabilities / total_prob

            self.position = np.random.choice([self.position, next_position], p=normalized_probabilities)

            # Update the position directly in the list
            existing_polymerase_positions[existing_polymerase_positions.index(previous_position)] = self.position

        if self.position >= end_of_replication_position:
            # Assuming delete returns the updated list
            existing_polymerase_positions = self.delete(self.position, existing_polymerase_positions)

        return existing_polymerase_positions

    def change_histones(self, chromatine,CenH_positions):
        if self.position not in CenH_positions:
            if 0 <= self.position < len(chromatine.histones) and chromatine.histones[self.position] == 'U' and np.random.random() < 0.5:
                chromatine.histones[self.position] = 'A'
            elif 0 <= self.position < len(chromatine.histones) and chromatine.histones[self.position] == 'M':
                chromatine.histones[self.position] = 'U'

# ---------------------------------------------------------------------------------- #

#                       Initialize chromatine and polymerases 

# ---------------------------------------------------------------------------------- #

chromatine = Chromatine(chromatine_size,CenH_positions)
polymerases = [Polymerase(chromatine) for _ in range(polymerase_count)]

# Track existing polymerase positions using a list to avoid duplicates
existing_polymerase_positions = [polymerase.position for polymerase in polymerases]

# Create an empty dataframe to store the counting lists
columns = ['Time Steps', 'Polymerase Count', 'Burst Size','Burst Frequency','Poly Proba movement', 'Count A']
# result_df = pd.DataFrame(columns=columns)

file = open('ScanCenHBURSTorNOBURSTfullcsvWriting.csv','w')
file.write('Time Steps, Polymerase Count, Burst Size, Burst Frequency, Poly Proba movement, Count A\n')


for num_poly_burst in np.arange(1,7,1):

    for burst_frequency in np.arange(0,1,0.05):

        for right_movement_probability in np.arange(1e-4,1.05,5e-2):

            print(f'Size burst : {num_poly_burst}, burst frequency : {burst_frequency}, RMP : {right_movement_probability}')

            # Polymerase movement probabilities
            left_movement_probability = 1/2



            # ---------------------------------------------------------------------------------- #

            #                               Simulation loop 

            # ---------------------------------------------------------------------------------- #

            for frame in range(simulation_steps):

                polymerase_positions = []  # Clear the polymerase_positions list
                polymerase_count = 0
                noisy_changes_count = 0
                enzyme_changes_count = 0

                # Adding polymerases with a burst
                if frame % int((1 - burst_frequency) * simulation_steps) == 0:
                    existing_polymerase_positions = chromatine.BurstPoly(adding_position,num_poly_burst,existing_polymerase_positions)
                
                polymerases = [Polymerase(chromatine, position=pos) for pos in existing_polymerase_positions]

                for polymerase in reversed(polymerases):
                    existing_polymerase_positions =  polymerase.move(chromatine,existing_polymerase_positions)
                    polymerase.change_histones(chromatine, CenH_positions)
                    polymerase_positions.append(polymerase.position)  # Append the current position


                # Change the next histones based on the influence of first neighbors
                position = np.random.randint(1, chromatine_size)
                
                # CenH_positions = chromatine.CenHRegion(CenH_positions,CenHSart,McenHDensity=MCenHDensity)
                # Use p_recruitment and p_change probabilities with decreasing probability with vicinity
                if np.random.random() < alpha:
                    enzyme_changes_count = chromatine.change_next_histones(position,CenH_positions, p_recruitment=recruitment_probability,
                                                                    p_change=change_probability, enzyme_changes=enzyme_changes_count,
                                                                    nth_neighbor=np.random.randint(1, chromatine_size))
                else:
                    noisy_changes_count = chromatine.noisy_transition(position,CenH_positions, noisy_transition_probability, noisy_changes_count)


                # Randomly add new polymerase at the beginning of the chromatine with a certain probability
                # if np.random.random() < chromatine.adding_poly_proba(adding_position, existing_polymerase_positions):
                    # Add new polymerases with non-overlapping random positions
                    # previous_poly_positions = polymerase_positions
                    # chromatine.add_polymerases(1,  adding_position, existing_polymerase_positions)
                    # new_polymerase_positions = existing_polymerase_positions[-1:]
                    # new_polymerases = [Polymerase(chromatine, position=pos) for pos in new_polymerase_positions]
                    # if previous_poly_positions != existing_polymerase_positions:
                        # polymerases.extend(new_polymerases)



                if frame%100 == 0:
                    # print(frame)
                    # Update the number of polymerases and active histones lists
                    polymerase_count = len(polymerases)
                    active_histone_count = np.sum(np.isin(chromatine.histones, ['M', 'A']))
                    acetylated_histone_count = np.sum(chromatine.histones == 'A')
                    methylated_histone_count = np.sum(chromatine.histones == 'M')
                    unmodified_histone_count = np.sum(chromatine.histones == 'U')
                    
                    chromatine_array = chromatine.ChromatineVisualisation()

                    count_A = np.count_nonzero(np.fromiter((nucleo == 'A' for nucleo in chromatine.histones[gene_position]), dtype=bool))


                    # Append data to the dataframe
                    file.write(f'{frame + 1}, {polymerase_count}, {num_poly_burst},{burst_frequency}, {right_movement_probability}, {count_A} \n')
                    
                    # result_df = pd.concat([result_df, pd.DataFrame([{'Time Steps': frame + 1,
                                                                # 'Polymerase Count': polymerase_count,
                                                                # 'Burst Size' : num_poly_burst,
                                                                # 'Burst Frequency': burst_frequency,
                                                                # 'Poly Proba movement' : right_movement_probability, 
                                                                # 'A in gene': count_A}])],ignore_index=True)
                    
print('Done')

file.close()

# Save the dataframe to a CSV file
   

current_directory = os.getcwd()

os.makedirs(current_directory, exist_ok=True)

csv_filename = os.path.join(current_directory, f'ScanCenHBURSTorNOBURST.csv')
print(csv_filename)
# result_df.to_csv(csv_filename, index=False)



# burst or non burst with the same number of polymerases during the sim