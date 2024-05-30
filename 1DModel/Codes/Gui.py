import tkinter as tk
from tkinter import ttk
import numpy as np
import pandas as pd

# Paste your functions and classes here

# ---------------------------------------------------------------------------------- #

#                                        Classes 

# ---------------------------------------------------------------------------------- #
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

    def noisy_transition(self, position,CenH_positions, noisy_transition_probability, noisy_changes, time, nucleosome_changes):
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

        return noisy_changes, time

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

    def change_next_histones(self,position ,CenH_positions, time, nucleosome_changes, p_recruitment, p_change, enzyme_changes, nth_neighbor):
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


        return enzyme_changes, time

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
    
    def CellCycle(self):
        num_unmodified = len(self.histones) // 2
        unmodified_positions = np.random.choice(len(self.histones), size=num_unmodified, replace=False)
        self.histones[unmodified_positions] = 'U'

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


class SimulationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Chromatine Simulation")

        # Create parameter input fields
        self.create_parameter_inputs()

        # Create start button
        start_button = ttk.Button(root, text="Start Simulation", command=self.start_simulation)
        start_button.pack(pady=10)

        # Output text area
        self.output_text = tk.Text(root, height=10, width=50)
        self.output_text.pack(pady=10)

    def create_parameter_inputs(self):
        parameters = [
            ("Chromatine Size", 198),
            ("Polymerase Count", 0),
            ("Simulation Steps", 30000),
            # Add more parameters here
        ]

        self.parameter_entries = {}
        for param_name, default_value in parameters:
            label = ttk.Label(self.root, text=param_name)
            label.pack(pady=5)
            entry = ttk.Entry(self.root)
            entry.insert(0, str(default_value))
            entry.pack(pady=5)
            self.parameter_entries[param_name] = entry

    def start_simulation(self):
        # Retrieve parameter values
        params = {param_name: int(entry.get()) for param_name, entry in self.parameter_entries.items()}

        # Run simulation with provided parameters
        result_df = self.run_simulation(**params)

        # Display results
        self.output_text.delete("1.0", tk.END)
        self.output_text.insert(tk.END, result_df)

    def run_simulation(self, chromatine_size, polymerase_count, simulation_steps):
        # Paste your simulation code here

        # For demonstration, let's return a DataFrame with some sample data
        result_df = pd.DataFrame({
            "Time Steps": np.arange(1, simulation_steps + 1),
            "Polymerase Count": np.random.randint(0, 10, size=simulation_steps),
            # Add more columns as needed
        })

        return result_df


if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationApp(root)
    root.mainloop()
