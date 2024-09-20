import re
import os
import numpy as np
from itertools import combinations

# dictionary for valence electrons
valence_electrons = {
    'H': 1, 'He': 2, 'O': 6, 'C': 4, 'N': 5, 'Na': 1, 'Cl': 7, 'S': 6, 'P': 5,
    # add more elements as needed
}

#generates list of  numbers
def generate_nums(max):
    return [num for num in range(0, max)]

#generates all unique combinations
def generate_pattern(numbers):
    return list(combinations(numbers, 2))

# asks for input
inp = input("Molecule (eg.H2O, NaCl): ")
type = input("Is this molecule ionic or covalent?: ").lower()

#seperates molecule into its individual elements
split = []
split = re.findall(r'([A-Z][a-z]*)(\d*)', inp)
newSplit = []
print(split)
for element, count in split:
    if count == '':
        count = 1
    else:
        count = int(count)
    newSplit.extend([element] * count)
print(newSplit)

#rearranges so central is first
central = min(newSplit, key=newSplit.count)
newSplit.insert(0, newSplit.pop(newSplit.index(central)))

#gets valence of each element
valence = []
for element in newSplit:
    if element in valence_electrons:
        valence.append(valence_electrons[element])
    else:
        raise ValueError(f"Unknown element {element}. Please add it to the valence_electrons dictionary.")
print(valence)

#creates a file
with open("model_temp.wcsp", "w+") as f:
    #write the header
    length = len(newSplit)
    nums = generate_nums(len(newSplit))
    combos = generate_pattern(nums)
    num_vars = length + len(combos) 
    f.write(f"Model {num_vars} 9 x 10\n")

    #writes variable domains
    domain = ""
    for i in valence:
        f.write(f"{i} ")
        domain += f"{i}"
    for x in combos:
        f.write("7 ")
        domain += "7"
    f.write("\n")

    #constraint for 2*bonds is var
    counter = 0
    index = domain.find("7")
    for i in range (index, num_vars):
        f.write(f"1 {i} 0 3\n1 10\n3 10\n5 10\n")
        counter += 1

    #add constraints for bonds according to molecule type
    if type == "covalent":
        matrix = np.zeros([length, length], dtype=int)
        for i in range(length):
            for j in range(length):
                if i == j:
                    matrix[i, j] = i
                else:
                    matrix[i, j] = i + (len(newSplit) - 1) + j
        for x in range(len(newSplit)):
            f.write(f"3  {' '.join(map(str, matrix[x]))} -1 wsum hard 10 == ")
            #octet and duet rule
            if newSplit[x] == 'H' or newSplit[i] == "He":
                f.write("2\n")
            else:
                f.write("8\n")
            counter += 1
    
    #add ionic

#reads file to replace placeholder
with open("model_temp.wcsp", "r") as f:
    file_contents = f.read()
file_contents = file_contents.replace('x', str(counter))

#creates new file with updated data
with open("model.wcsp", "w") as f:
    f.write(file_contents)
os.remove("model_temp.wcsp")
print("WCSP model written to 'model.wcsp'.")