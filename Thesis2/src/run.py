import subprocess

# Replace 'your_program' with the actual name of your compiled C program
program_name = './a.out'

# Add any command-line arguments your C program needs
arguments = ['ATTGACAACCGCGTCCGCCGCCGCTTCAAGGGCCAGTACTTGATGCCCAACATTGGCTACGGCTCCAACAAGCGCACCCGCCACATGTTGCCCACCGGCT', 'IDNRVRRRFKGQYLMPNIGYGSNKRTRHMLPTGF']

# Run the C program using subprocess
process = subprocess.Popen([program_name] + arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# Wait for the process to finish and capture the output
stdout, stderr = process.communicate()

# Print the output
print(stdout.decode())
# print("Output:", stdout.decode())
# print("Error:", stderr.decode())