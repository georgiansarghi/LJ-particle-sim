def read_xyz_file(file_name):
    """Reads an xyz file and returns a dictionary with time as keys and particle data as values."""
    with open(file_name, 'r') as file:
        lines = file.readlines()

    data = {}
    i = 0
    while i < len(lines):
        particle_count = int(lines[i].strip())
        time_stamp = lines[i + 1].strip()
        particles = lines[i + 2:i + 2 + particle_count]
        data[time_stamp] = [p.strip() for p in particles]
        i += 2 + particle_count

    return data

def merge_xyz_files(file_names):
    """Merges multiple xyz files."""
    all_data = {}

    # Read and merge data from each file
    for file_name in file_names:
        data = read_xyz_file(file_name)
        for time_stamp, particles in data.items():
            if time_stamp not in all_data:
                all_data[time_stamp] = particles
            else:
                all_data[time_stamp].extend(particles)

    return all_data

def write_merged_file(merged_data, output_file_name):
    """Writes the merged data to a new file."""
    with open(output_file_name, 'w') as file:
        for time_stamp, particles in merged_data.items():
            file.write(str(len(particles)) + '\n')
            file.write(time_stamp + '\n')
            file.write('\n'.join(particles) + '\n')

# Usage
file_names = ['S5.0000.xyz', 'S5.0001.xyz', 'S5.0002.xyz', 'S5.0003.xyz']
merged_data = merge_xyz_files(file_names)
write_merged_file(merged_data, 'merged_output.xyz')
