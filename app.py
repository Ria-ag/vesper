import os
import re
import subprocess
import numpy as np
from flask_cors import CORS
from itertools import combinations
from flask import Flask, request, render_template, jsonify
from chemlib import Element

app = Flask(__name__)
CORS(app)

def createWCSPFile(inp, btype):
    # dictionary for valence electrons
    valence_electrons = {
        'H': 1, 'He': 2, 'O': 6, 'C': 4, 'N': 5, 'Na': 1, 'Cl': 7, 'S': 6, 'P': 5, 'B': 3, 'F': 7,
        # replace this eventually with pulling this data from library
    }
    electronegativities = {
        'H': 2.20, 'He': None, 'O': 3.44, 'C': 2.55, 'N': 3.04, 'Na': 0.93, 'Cl': 3.16, 'S': 2.58, 'P': 2.19, 'B': 2.04, 'F': 3.98,
    }

    #generates list of  numbers
    def generate_nums(max):
        return [num for num in range(0, max)]

    #generates all unique combinations
    def generate_pattern(numbers):
        return list(combinations(numbers, 2))

    #seperates molecule into its individual elements
    split = []
    split = re.findall(r'([A-Z][a-z]*)(\d*)', inp)
    newSplit = []
    for element, count in split:
        if count == '':
            count = 1
        else:
            count = int(count)
        newSplit.extend([element] * count)

    #gets electronegativity of each element
    electronegativity = []
    for element in newSplit:
        if element in electronegativities:
            electronegativity.append(electronegativities[element])
        else:
            raise ValueError(f"Unknown element {element}. Please add it to the electronegativities dictionary.") #will not need once replaced with library
        
    #find central atoms
    excluded = {"H", "F", "Cl", "Br", "I"}
    indexx = 0
    finalIndex = 0
    potentialCenrals = []
    minElec = float("inf")
    for indexx, element in enumerate(newSplit):
        if element not in excluded and electronegativity[indexx] < minElec:
            minElec = electronegativity[indexx]
            potentialCentral = [indexx]
            finalIndex = indexx
        indexx += 1

    central = []
    noncentral = []
    for element in newSplit:
        if element == newSplit[finalIndex]:
            central.append(element)
        else:
            noncentral.append(element)

    newSplit = central + noncentral

     #gets valence of each element
    totalVal = 0
    valence = []
    for element in newSplit:
        if element in valence_electrons:
            valence.append(valence_electrons[element])
            totalVal += valence_electrons[element]
        else:
            raise ValueError(f"Unknown element {element}. Please add it to the valence_electrons dictionary.") #will not need once replaced with library

    #creates a file
    with open("model_temp.wcsp", "w+") as f:
        #write the header
        length = len(newSplit)
        nums = generate_nums(len(newSplit))
        combos = generate_pattern(nums)
        num_vars = length + len(combos) 
        f.write(f"Model {num_vars} 9 x 10\n")
        btype = btype.lower()

        #writes variable domains
        domain = ""
        for i in valence:
            if btype == "ionic" and i > 4:
                f.write("9 ")
                domain += "9"
            else:
                f.write(f"{i} ")
                domain += f"{i}"
        for x in combos:
            f.write("7 ")
            domain += "7"
        f.write("\n")
    
        #create matrix for each element pairs
        matrix = np.zeros([length, length], dtype=int)
        next = 1
        for i in range(length):
            for j in range(length):
                if i == j:
                    matrix[i, j] = i
                else:
                    if i == 0:
                        matrix[i, j] = i + (length - 1) + j
                    else:
                        if matrix.item((j, i)) != 0:
                            matrix[i, j] = matrix.item((j, i))
                        else: 
                            matrix[i, j] = matrix.item((0, length-1)) + next
                            next += 1

        #add constraints for bonds according to molecule type
        counter = 0
        if btype == "covalent":
             #constraint for 2*bonds is var
            for i in range (length, num_vars):
                f.write(f"1 {i} 0 3\n1 10\n3 10\n5 10\n")
                counter += 1
            for x in range(len(newSplit)):
                f.write(f"{length} {' '.join(map(str, matrix[x]))} -1 wsum hard 10 == ")
                #octet and duet rule
                if newSplit[x] == 'H' or newSplit[x] == "He":
                    f.write("2\n")
                else:
                    f.write("8\n")
                counter += 1
            #has to equal total valence count
            f.write(f"{num_vars} ")
            for i in range (num_vars):
                f.write(f"{i} ")
            f.write(f"-1 wsum hard 10 == {totalVal}\n")
            counter += 1
            #central atom chain
            indexxx = length
            k = len(central)
            if k > 1:
               for i in range(len(central) - 1):
                    f.write(f"1 {indexxx} 0 1\n0 10\n")
                    indexxx += k
                    k -= 1
                    counter += 1
            #each must bond to at least one central atom
            index = length + len(central) - 1
            if len(central) > 1:
                for i in range(len(noncentral)):
                    f.write(f"{len(central)} ")
                    for j in range(len(central)):
                        constraintIndex = (index + (len(noncentral) * j))
                        if (j > 0 and len(central) > 2):
                            constraintIndex += 1
                        f.write(f"{constraintIndex} ")
                    f.write(f"-1 wsum hard 10 > 0\n")
                    index += 1
                    counter += 1
            #constraint for lone pairs to harm the molecule
            for i in range (length):
                f.write(f"1 {i} 10 1\n0 0\n")
            #Either bonded to O or a central atom
            if k == 1:
                HIndex = 0
                OIndex = 0
                for i in newSplit:
                    if i == "O":
                        OIndex = newSplit.index(i)
                        if OIndex != 0:
                            for i in newSplit:
                                if i == "H":
                                    HIndex = newSplit.index(i)
                                    f.write(f"2 {matrix[HIndex][OIndex]} {matrix[OIndex][0]} -1 wsum hard 10 != 0\n")
                                    counter += 1

        else:
            for i in range (length, num_vars):
                f.write(f"1 {i} 10 1\n0 0\n")
                counter += 1
            i = 0
            for x in valence:
                if x >= 4:
                    f.write(f"{length} {' '.join(map(str, matrix[i]))} -1 wsum hard 10 == 8\n")
                else:
                    #looses electrons
                    f.write(f"{length} {' '.join(map(str, matrix[i]))} -1 wsum hard 10 == 0\n")
                i += 1
                counter += 1

    #reads file to replace placeholder
    with open("model_temp.wcsp", "r") as f:
        file_contents = f.read()
    file_contents = file_contents.replace('x', str(counter))

    #creates new file with updated data
    with open("model.wcsp", "w") as f:
        f.write(file_contents)
    os.remove("model_temp.wcsp")
    print("WCSP model written to 'model.wcsp'.")

    centrals = len(central)

    return "model.wcsp", num_vars, newSplit, centrals

#route to render the html page
@app.route('/')
def index():
    return render_template('Vesper.html')

#route to run the wcsp code
@app.route('/run_code', methods=['POST'])
def run_code():
    inp = request.form['inp']
    btype = request.form['type']
    
    #run og python code to generate wcsp file
    wcspFilePath, num_vars, newSplit, centrals = createWCSPFile(inp, btype)

    #ssh into the instance and run toulbar2
    ssh_command = f"toulbar2 model.wcsp -s -a"
    ssh_result = subprocess.run(ssh_command, shell=True, capture_output=True, text=True)

    #parse output to get all solutions
    output = ssh_result.stdout.splitlines()
    solutions = []
    for line in output:
        if "solution" in line:
            solutions.append(line.strip())

    #return or error
    if solutions:
        return jsonify({'output': solutions, 'num_vars': num_vars, 'elements': newSplit, 'centrals': centrals})
    else:
        return jsonify({'output': 'Solution not found'}), 500
    
if __name__ == '__main__':
    app.run(host='0.0.0.0')