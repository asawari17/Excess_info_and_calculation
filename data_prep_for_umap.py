import os
import pandas as pd

# Path to the folder containing all the folders with text files
main_folder = '/Users/asawaripagare/Desktop/all_data/'

# Path to the output Excel file
output_file = '/Users/asawaripagare/Desktop/i_features.xlsx'

# Define a mapping of R1_1 values to numbers
r1_1_mapping = {
    'R1_1': 10, 'R1_2': 1, 'R1_3': 12, 'R1_4': 11, 'R1_6': 2, 'R1_7': 3, 'R1_9': 13,
    'R2a_1': 14, 'R2a_2': 17, 'R2a_3': 15, 'R2b_1': 18, 'R2b_2': 21, 'R2b_3': 20, 'R2c_1': 4, 'R2c_2': 5,
    'R2c_4': 6, 'R2c_6': 16, 'R3a_1': 7, 'R3a_4': 8,
    'R3a_7': 9, 'R3b_1': 23, 'R3b_4': 22,
    'R3b_7': 24, 'R3c_1': 19
}

# Get the folder names from the mapping
folders_to_process = set(r1_1_mapping.keys())

# Initialize an empty list to store the data
data = []

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
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                    col1, col2 = [], []
                    for line in lines:
                        try:
                            values = line.split()
                            col1.append(float(values[0]))
                            col2.append(float(values[1]))
                        except (ValueError, IndexError):
                            # Skip lines that do not contain valid numeric data
                            continue
                    if not col1 or not col2:
                        continue  # Skip if col1 or col2 is empty
                    argmax_col2_over_max_col1 = col2.index(max(col2)) / col1.index(max(col1))
                    max_col2 = max(col2)
                    col2_series = pd.Series(col2)  # Convert col2 to a pandas Series
                    first_entry_over_max_entry = col2_series.iloc[0] / max_col2
                    last_entry_over_max_entry = col2_series.iloc[-1] / max_col2
                # Extract information from file name
                file_name_parts = file_name.split('_')
                modified_index = int(file_name_parts[file_name_parts.index('modified') + 1])
                r1_1_value = r1_1_mapping.get(folder_name, None)
                if r1_1_value is None:
                    continue  # Skip if r1_1_value is None
                # Append the values for this file to the data list
                data.append([argmax_col2_over_max_col1, max_col2, first_entry_over_max_entry, last_entry_over_max_entry, modified_index, r1_1_value])

# Create a DataFrame from the data
df = pd.DataFrame(data, columns=['attr1', 'attr2', 'attr3', 'attr4', 'itr_number', 'Node_identity'])

# Write the DataFrame to an Excel file
df.to_excel(output_file, index=False)
