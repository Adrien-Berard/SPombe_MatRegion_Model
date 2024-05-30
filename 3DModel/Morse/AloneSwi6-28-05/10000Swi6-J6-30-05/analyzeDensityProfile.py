from mpi4py import MPI
import pandas as pd

def read_profile(file_path):
    # Read the file
    with open(file_path, 'r') as file:
        # Read the metadata lines
        meta_lines = []
        while True:
            line = file.readline()
            if not line.startswith('#'):
                meta_lines.append(line)
                break
            meta_lines.append(line)
        
        # Read the rest of the file into a pandas DataFrame
        data = pd.read_csv(file_path, delim_whitespace=True, comment='#', header=None, skiprows=len(meta_lines))

    # Extract metadata
    timestep = int(meta_lines[-1].split()[0])
    number_of_chunks = int(meta_lines[-1].split()[1])
    total_count = int(meta_lines[-1].split()[2])

    # Assign column names to the DataFrame based on the provided format
    data.columns = ['Chunk', 'Coord1', 'Coord2', 'Coord3', 'Ncount', 'Density']

    return timestep, number_of_chunks, total_count, data

def process_data(data, rank, size):
    # Split the data based on rank
    chunk_size = len(data) // size
    start = rank * chunk_size
    end = (rank + 1) * chunk_size if rank != size - 1 else len(data)
    local_data = data[start:end]

    # Perform some processing, e.g., calculating the average density
    avg_density = local_data['Density'].mean()

    return avg_density

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    file_path = 'density3d.profile'
    if rank == 0:
        timestep, number_of_chunks, total_count, data = read_profile(file_path)
    else:
        data = None

    # Broadcast data to all processes
    data = comm.bcast(data, root=0)

    # Each process handles its portion of the data
    local_avg_density = process_data(data, rank, size)

    # Gather results from all processes
    all_avg_densities = comm.gather(local_avg_density, root=0)

    if rank == 0:
        global_avg_density = sum(all_avg_densities) / size
        print(f"Global average density: {global_avg_density}")

if __name__ == "__main__":
    main()
