#截取有用部分
#K00097_**_eco:b0052     GCF_000190995.1_ASM19099v1_**_WP_000241242.1
#GCF_000190995.1         K00097
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input file content
with open(input_file, "r") as file:
    lines = file.readlines()

new_lines = []

# Process each line
for line in lines:
    parts = line.strip().split('\t')
    if len(parts) >= 2:
        col1, col2 = parts[0], parts[1]
        col1 = col1.split("_**_")[0]        
        # Remove content after the second _ in the second column
        col2 = "_".join(col2.split("_", 2)[:2])
        
        # Swap positions of the first and second columns
        new_line = "{}\t{}\t{}\n".format(col2, col1, '\t'.join(parts[2:]))
        new_lines.append(new_line)

# Write the processed content to the output file
with open(output_file, "w") as file:
    file.writelines(new_lines)

print("Processing complete. Results have been written to", output_file)

