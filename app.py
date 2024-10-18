import os
import re
import subprocess
import numpy as np
from flask_cors import CORS
from itertools import combinations
from flask import Flask, request, render_template, jsonify

app = Flask(__name__)
CORS(app)

def createWCSPFile(inp, type):
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

    #creates a file
    with open("model_temp.wcsp", "w+") as f:
        #write the header
        length = len(newSplit)
        nums = generate_nums(len(newSplit))
        combos = generate_pattern(nums)
        num_vars = length + len(combos) 
        f.write(f"Model {num_vars} 9 x 10\n")
        type = type.lower()

        #writes variable domains
        domain = ""
        for i in valence:
            if type == "ionic" and i > 4:
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
        if type == "covalent":
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
            index = length
            while index < length + (length - 1):
                f.write(f"1 {index} 0 1\n0 10\n")
                counter += 1
                index += 1
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

    return "model.wcsp", num_vars, newSplit

#route to render the html page
@app.route('/')
def index():
    return render_template('Vesper.html')

#route to run the wcsp code
@app.route('/run_code', methods=['POST'])
def run_code():
    inp = request.form['inp']
    type = request.form['type']
    
    #run og python code to generate wcsp file
    wcspFilePath, num_vars, newSplit = createWCSPFile(inp, type)

    #ssh into the instance and run toulbar2
    ssh_command = f"toulbar2 model.wcsp -s -a"
    ssh_result = subprocess.run(ssh_command, shell=True, capture_output=True, text=True)

    #parse output to get solution
    output = ssh_result.stdout.splitlines()
    solution = None
    for line in output:
        if "solution" in line:
            solution = line.strip()
            break

    #return or error
    if solution:
        return jsonify({'output': solution, 'num_vars': num_vars, 'elements': newSplit})
    else:
        return jsonify({'output': 'Solution not found'}), 500
    
if __name__ == '__main__':
    app.run(host='0.0.0.0')