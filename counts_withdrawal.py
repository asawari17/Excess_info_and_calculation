import os
import pandas as pd

# Path to the folder containing all the folders with text files
main_folder = 'path/to/folder'

# Path to the output Excel file
output_counts_file = 'output_file.xlsx'

# Define a mapping of R1_1 values to numbers
mapping = {
    'R1_1': 10, 'R1_2': 1, 'R1_3': 12, 'R1_4': 11, 'R1_6': 2, 'R1_7': 3, 'R1_9': 13,
    'R2a_1': 14, 'R2a_2': 17, 'R2a_3': 15, 'R2b_1': 18, 'R2b_2': 21, 'R2b_3': 20, 'R2c_1': 4, 'R2c_2': 5,
    'R2c_4': 6, 'R2c_6': 16, 'R3a_1': 7, 'R3a_4': 8,
    'R3a_7': 9, 'R3b_1': 23, 'R3b_4': 22,
    'R3b_7': 24, 'R3c_1': 19
}

# Get the folder names from the mapping
folders_to_process = set(mapping.keys())

# Initialize an empty list to store the data
data = []

# Initialize a dictionary to count the occurrences where attr1 != 1
counts = {r1_1_value: 0 for r1_1_value in mapping.values()}

# Iterate through each folder in the main folder
for folder_name in os.listdir(main_folder):
    if folder_name not in folders_to_process:
        continue
    folder_path = os.path.join(main_folder, folder_name)
    if os.path.isdir(folder_path):
        # Iterate through each file in the folder
        for file_name in os.listdir(folder_path):
            if file_name.endswith('.txt'):
                file_path = os.path.join(folder_path, file_name)
                argmax_col2_over_max_col1 = None
                max_col2_over_last_entry = None
                max_minus_last = None
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                    col1 = [float(line.split()[0]) for line in lines]
                    col2 = [float(line.split()[1]) for line in lines]
                    argmax_col2_over_max_col1 = col2.index(max(col2)) / (col1.index(max(col1)))
                # Extract information from file name
                file_name_parts = file_name.split('_')
                modified_index = int(file_name_parts[file_name_parts.index('modified') + 1])
                r1_1_value = None
                if len(file_name_parts) >= 3:
                    r1_1_part = file_name_parts[-3] + '_' + file_name_parts[-2]
                    if r1_1_part in mapping:
                        r1_1_value = mapping[r1_1_part]

                # Update the count if attr1 != 1
                if argmax_col2_over_max_col1 != 1 and r1_1_value is not None:
                    counts[r1_1_value] += 1


# Create a DataFrame from the counts dictionary
counts_df = pd.DataFrame(list(counts.items()), columns=['Node_identity', 'attr1_not_1_count'])

# Write the counts DataFrame to an Excel file
counts_df.to_excel(output_counts_file, index=False)
