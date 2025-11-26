import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc
from math import gcd
from functools import reduce
# from sklearn.model_selection import train_test_split
# from sklearn.tree import DecisionTreeRegressor
# from sklearn import tree

def generate_6_three_digit_ternary():
    ternary_numbers = []
    for i in range(6):
        ternary_number=[0,0,0]
        val = (-1)**(i%2+1)
        ternary_number[int(i/2)] = val
        ternary_numbers.append(ternary_number)
    return ternary_numbers

def generate_all_three_digit_ternary():
    ternary_numbers = []
    for i in range(3**3-1,-1,-1):
        resid = i
        ternary_number = [0,0,0]
        for j in range(3):
            power = 2-j
            ternary = int(resid/(3**power))
            resid = resid - ternary*(3**power)
            ternary_number[j] = ternary

        ternary_numbers.append(ternary_number)

    return ternary_numbers

def calculate_common_divisors(numbers):

    def find_divisors(n):
        divisors = []
        for i in range(1, int(n**0.5) + 1):
            if n % i == 0:
                divisors.append(i)
                if i != n // i:
                    divisors.append(n // i)
        return set(divisors)

    gcd_value = reduce(gcd, numbers)
    return sorted(find_divisors(gcd_value))

def calculate_strain_ratios(digit_numbers):
    strain_ratios = []
    maxv = max(max(digit_numbers))
    for digit_number in digit_numbers:
        vect = np.array(digit_number) - maxv/2.0
        norm = np.linalg.norm(vect)
        if norm > 0:
            strain_ratio = vect / norm
            strain_ratios.append(strain_ratio)
#   with open("./results/strain_ratios.txt", 'w') as f:
#       for strain_ratio in strain_ratios:
#           f.write(f"{strain_ratio}\n")
    return strain_ratios


def latin_hypercube_sampling(num_data, dimension):
    sampler = qmc.LatinHypercube(d=dimension)
    sample = sampler.random(n=num_data)
    return sample


def process_lhs_samples(lhs_samples, strain_ratios):
    for lhs_sample in lhs_samples:
        strain_ratio = 2.0 * np.array(lhs_sample) - 1.0
        norm = np.linalg.norm(strain_ratio)
        if norm > 0:
            strain_ratio = strain_ratio / norm
            strain_ratios.append(strain_ratio)

    # Print the strain ratios and check normalization
    for ipath, strain_ratio in enumerate(strain_ratios):
        if np.linalg.norm(strain_ratio) - 1.0 > 1e-6:
            print("Error in normalization")
            break

def plot_strain_ratios(strain_ratios,filename):
    # Generate 30 graphs by plotting each pair of columns
    num_columns = strain_ratios.shape[1]
    graph_count = 0
    plt.figure(figsize=(9, 9))
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    for i in range(num_columns):
        for j in range(num_columns):
            graph_count += 1
            plt.subplot(3, 3, graph_count)
            plt.scatter(strain_ratios[:, i], strain_ratios[:, j], s=1, alpha=0.5)
            plt.xlabel(f"Col {i}")
            plt.ylabel(f"Col {j}")
            plt.title(f"{i} vs {j}")
            plt.grid(True)
    plt.savefig(filename, bbox_inches='tight', dpi=300)

def generate_6_strain_ratios():
    ternary_numbers = generate_6_three_digit_ternary()
    strain_ratios = np.array(ternary_numbers)

    return strain_ratios

def generate_and_plot_strain_ratios(npath_total):
    # Generate and print a three-digit ternary number
    digits = generate_all_three_digit_ternary()
    # Make uniformly distributed 3D strain ratios
    strain_ratios = calculate_strain_ratios(digits)
    # Make randomly distributed 3D strain ratios
    npath_random = npath_total - len(strain_ratios)
    lhs_samples = latin_hypercube_sampling(npath_random, 3)
    # Call the function
    process_lhs_samples(lhs_samples, strain_ratios)
    # Convert strain_ratios to a numpy array for easier manipulation
    strain_ratios = np.array(strain_ratios)
    # Call the function
    # plot_strain_ratios(strain_ratios,"./figures/scatter_plot_parametric.png")

    return strain_ratios


def convert_to_physical_strains(strain_ratios_parametric, max_strains, min_strains):
    strain_ratios = np.array(strain_ratios_parametric, dtype=float)
    for ipath in range(len(strain_ratios_parametric)):
        for icomp in range(3):
            if abs(strain_ratios_parametric[ipath][icomp]) < 1.0e-6:
                strain_value = 0.0
            elif strain_ratios_parametric[ipath][icomp] > 0.0:
                strain_value = max_strains[icomp] * strain_ratios_parametric[ipath][icomp]
            elif strain_ratios_parametric[ipath][icomp] < 0.0:
                strain_value = -min_strains[icomp] * strain_ratios_parametric[ipath][icomp]
            strain_ratios[ipath][icomp] = strain_value

    # plot_strain_ratios(strain_ratios,"./figures/scatter_plot_physical.png")
    return strain_ratios

def strain_matrix(voigt, idm):
    strain = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            strain[i, j] = voigt[idm[i][j]-1]
            if i != j:
                strain[i, j] = strain[i, j] / 2.0
    return strain

def matrix_to_voigt(mat):
    voigt = np.zeros(6)
    idm = [[1, 4, 6],
           [4, 2, 5],
           [6, 5, 3]]
    for i in range(3):
        for j in range(3):
            voigt[idm[i][j]-1] = mat[i, j]
    return voigt

def compute_stress_strain(npath_total, nstep, strain_ratios, yng, poi):
    idm = [[1, 4, 6],
           [4, 2, 5],
           [6, 5, 3]]
    eye = np.zeros((3, 3))
    for i in range(3):
        eye[i, i] = 1.0

    bulk = yng / (3.0 * (1.0 - 2.0 * poi))
    shear = yng / (2.0 * (1.0 + poi))
    lame1 = bulk - (2.0 / 3.0) * shear
    lame2 = shear

    stress_strains = np.zeros((npath_total, nstep, 12))

    for ipath, strain_ratio in enumerate(strain_ratios):

        for istep in range(nstep):
            factor_inc = float(istep + 1)/float(nstep)
            strain_voigt = np.zeros(6)
            strain_voigt[0] = factor_inc * strain_ratio[0]
            strain_voigt[1] = factor_inc * strain_ratio[1]
            strain_voigt[3] = factor_inc * strain_ratio[2]
            strain = strain_matrix(strain_voigt, idm)
            stress = lame1 * np.trace(strain) * eye + 2.0 * lame2 * strain
            stress_voigt = matrix_to_voigt(stress)

            stress_strains[ipath, istep, :6] = strain_voigt
            stress_strains[ipath, istep, 6:] = stress_voigt

    return stress_strains

def create_dataset(stress_strains, npath_total, nstep):
    dataset = np.zeros((npath_total * nstep, 12))
    for ipath in range(npath_total):
        for istep in range(nstep):
            dataset[ipath * nstep + istep, :] = stress_strains[ipath, istep, :]
    return dataset

def plot_stress_strain_relationship(dataset, filename="./figures/stress_strain_plots.png"):
    figure = plt.figure(figsize=(9, 9))
    figure.subplots_adjust(hspace=0.5, wspace=0.5)
    comps = [0,1,3]
    for i,icomp in enumerate(comps):
        for j,jcomp in enumerate(comps):
            axis = figure.add_subplot(3, 3, i * 3 + j + 1)
            axis.plot(dataset[:, icomp], dataset[:, jcomp + 6], 'o', markersize=1, alpha=0.5)
            axis.set_xlabel(f"Strain {icomp+1}")
            axis.set_ylabel(f"Stress {jcomp+1}")
            axis.grid(True)
    plt.savefig(filename, bbox_inches='tight', dpi=300)

def save_and_plot_datasets(dataset, npath_total, nstep):
    with open("./results/dataset.txt", 'w') as f:
        for data in dataset:
            outline = f"{data[0]:>20.15e}  {data[1]:>20.15e}  {data[3]:>20.15e}  {data[6]:>20.15e}  {data[7]:>20.15e}  {data[9]:>20.15e}"
            f.write(outline + "\n")
    plot_stress_strain_relationship(dataset)

## Conditions
def read_simulation_parameters(input_file):
    with open(input_file, 'r') as file:
        params = {}
        for line in file:
            key, value = line.strip().split('=')
            key = key.strip()
            value = value.strip()
            if ',' in value:
                params[key] = [float(v) for v in value.split(',')]
            else:
                params[key] = float(value) if '.' in value or 'e' in value.lower() else int(value)
    return params

# Main part

## Read parameters from input file
input_file = "./condition_datamake.txt"
params = read_simulation_parameters(input_file)
npath_total = params['npath_total']
nstep = params['nstep']
max_strains = [ params['max_strain1'], params['max_strain2'], params['max_strain3'],
                params['max_strain4'], params['max_strain5'], params['max_strain6'] ]
min_strains = [ params['min_strain1'], params['min_strain2'], params['min_strain3'],
                params['min_strain4'], params['min_strain5'], params['min_strain6'] ]
yng = params['yng']
poi = params['poi']

## Generate dataset
if npath_total == 6:
    strain_ratios_parametric = generate_6_strain_ratios()
else:
    strain_ratios_parametric = generate_and_plot_strain_ratios(npath_total)

strain_ratios = convert_to_physical_strains(strain_ratios_parametric, max_strains, min_strains)
stress_strains = compute_stress_strain(npath_total, nstep, strain_ratios, yng, poi)
dataset = create_dataset(stress_strains, npath_total, nstep)
save_and_plot_datasets(dataset, npath_total, nstep)

## Decision tree
# generate_decision_trees_and_extract_leaf_data(dataset)
