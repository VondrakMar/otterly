from ase.io import read,write
from ase.visualize import view
import matplotlib.pyplot as plt
import numpy as np
from matscipy.neighbours import neighbour_list
from copy import deepcopy
from collections import defaultdict
import pandas as pd

exp_fuck = [0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,4.0,5.0]

def get_water_charges(mols,charge_keyword="kqeq_charges",ion_index=None):
    if ion_index is None:
        N_mols = len(mols[0])//3
        qs = {i: [] for i in range(N_mols)}
        for mol in mols:
            mol_qs = mol.arrays[charge_keyword]
            sum_qs = np.sum(mol_qs.reshape(-1,3),axis=1)
            for i in range(N_mols):
                qs[i].append(sum_qs[i])
    return qs

def plot_molecule_charges(mols,charge_keyword="kqeq_charges",ref_charges_keyword="dft_hirshfeld",ion_index = None,show_plot = True):
    '''
    This is the main function, here I can take list of ase.atoms, and if the ion_index is not None, it will create tmp mols where I remove ion from the structures. Then it return x, y, y_ref used for plotting.
    '''
    assert all(len(mol) == len(mols[0]) for mol in mols), "Your list is not of the same sized clusters, this code was written solely for charges of the scaled water"
    if ion_index is not None:
        tmp_mols = []
        ion_index.sort(reverse=True)
        for mol in mols:
            tmp_mol = deepcopy(mol)
            for ion in ion_index:
                del(tmp_mol[ion])
            tmp_mols.append(tmp_mol)
        qs = get_water_charges(tmp_mols,charge_keyword)
        qs_ref = get_water_charges(tmp_mols,ref_charges_keyword)
    elif ion_index is None:
        assert len(mols[0])//3 == len(mols[0])/3, "It seems you do not have correct number of atoms for pure water"
        tmp_mols = [deepcopy(mol) for mol in mols] # just so the variable name is the same as before
        qs = get_water_charges(tmp_mols,charge_keyword)
        qs_ref = get_water_charges(tmp_mols,ref_charges_keyword)
    N_mols = len(tmp_mols[0])//3
    for i in range(N_mols):
        plt.plot(exp_fuck,qs[i])
        plt.plot(exp_fuck,qs_ref[i],"--")
    if show_plot:
        plt.show()
    qs = np.array(qs).T
    print("qs",qs)
    qs_ref = np.array(qs_ref).T
    return exp_fuck,qs,qs_ref


def guess_ion(mol,rcut=1.5):
    tmp_mol = deepcopy(mol)
    tmp_mol.center(vacuum=10)
    i, j = neighbour_list('ij', tmp_mol, cutoff=rcut)
    syms = tmp_mol.symbols
    neighbors = defaultdict(list)
    for a, b in zip(i, j):
        neighbors[a].append(b)
        # neighbors[b].append(a)
    oxygen_with_wrong_hcount = []
    for idx, atom in enumerate(syms):
        if atom == 'O':
            h_count = sum(syms[n] == 'H' for n in neighbors[idx])
            if h_count != 2:
                oxygen_with_wrong_hcount.append(idx)
    assert len(oxygen_with_wrong_hcount) == 1, "this code is intended only for structures with one ion. Or the ion is not distinquidfgshable with your cuttoff"
    if len(oxygen_with_wrong_hcount) == 0:
        print("No ion has been found")
        return 0
    else:
        result = []
        result.append(oxygen_with_wrong_hcount[0])
        result.extend([i for i in neighbors[oxygen_with_wrong_hcount[0]]])
        return(result)


def save_csv(x,y,y_ref,file_name):
    '''
    This is a function only for saving stuff from the functions above
    '''
    df = pd.DataFrame({'exp_factor': x})
    print("hey",y_ref[0])
    for i in y_ref:
        df[f'hirsh_{i+1}'] = y_ref[i]
    df.to_csv(file_name, index=False)

        
# def save_csv(res_dict,file_name):
#     '''
#     This is a function only for saving stuff from the functions above
#     '''
#     df = pd.DataFrame({'exp_factor': res_dict["x"]})
#     kqeq_charges = []
#     for id_exp, exp_list in enumerate(res_dict["y"]):
#         tmp_charges = []
#         print(exp_list)
#         for id_mol, mol in enumerate(res_dict["y"][exp_list]):
#             print(mol)
#             tmp_charges.append(mol)
#         kqeq_charges.append(tmp_charges)
#     kqeq_charges = np.array(kqeq_charges).T
#     for i, kqeq in enumerate(kqeq_charges):
#         df[f'kqeq_{i+1}'] = kqeq
#     for i, hirsh in enumerate(res_dict["y_ref"]):
#         df[f'hirsh_{i+1}'] = hirsh    
#     df.to_csv(file_name, index=False)


neut = read("finalNEUT.xyz@:",format="extxyz")
#qs = get_water_charges(neut)
x,y,y_ref = plot_molecule_charges(neut,show_plot=False)
save_csv(x,y,y_ref,"neut_res.csv")


neg = read("finalNEG.xyz@:",format="extxyz")
neg_ion_id = guess_ion(neg[-1])
neg_res = plot_molecule_charges(neg,ion_index=neg_ion_id)
save_csv(neg_res,"neg_res.csv")

pos = read("finalPOS.xyz@:",format="extxyz")
pos_ion_id = guess_ion(pos[-1])
pos_res = plot_molecule_charges(pos,ion_index=pos_ion_id)
save_csv(pos_res,"pos_res.csv")




