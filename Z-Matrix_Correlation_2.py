import mdtraj as md
import sys
import numpy as np
import pandas as pd

## LOAD TRAJECTORY ##
OutFN = ["./ala1_DoFs.csv", "./ala1_covmat.csv"]
TrajIn = md.load("ala1/ala1_Prod_Restraint_Global.dcd", top="ala1/ala1.prmtop")
print (TrajIn)
bondgraph = TrajIn.topology.to_bondgraph()
bond_Adjacency = bondgraph.adj

## Number the atoms ##
SingleEdge = {}
MultipEdge = {}

counter = 0
for AtomIx in range(TrajIn.n_atoms):
    if (len(bond_Adjacency[TrajIn.topology.atom(AtomIx)]) > 1):
        MultipEdge[counter] = AtomIx
        counter += 1
        #print ("Bing!")
for AtomIx in range(TrajIn.n_atoms):
    if (len(bond_Adjacency[TrajIn.topology.atom(AtomIx)]) == 1):
        SingleEdge[counter] = AtomIx
        counter += 1
        #print ("Bong!")

AllEdges = {**MultipEdge, **SingleEdge}

MultipEdge_val = list(MultipEdge.values())
SingleEdge_val = list(SingleEdge.values())

##  Generate The Z-Matrix (empty) and fill the first 3 lines
zMatrix = pd.DataFrame(index=[x for x in range(TrajIn.n_atoms)], columns=["i","B","A","T"])
zMatrix.loc[0] = [1, None, None, None]
zMatrix.loc[1] = [2, 1, None, None]
zMatrix.loc[2] = [3, 2, 1, None]
##  Iterate over MultipEdge Atoms   ##
for AtomIx in range(3, len(MultipEdge)):
    print (TrajIn.topology.atom(MultipEdge[AtomIx]))
    zMatrix["i"][AtomIx] = AtomIx + 1
    for j in range(AtomIx-1, -1, -1):
        if (TrajIn.topology.atom(MultipEdge[j]) in bond_Adjacency[TrajIn.topology.atom(MultipEdge[AtomIx])]):
            zMatrix["B"][AtomIx] = j+1
            break
    for k in range(j-1, -1, -1):
        if (TrajIn.topology.atom(MultipEdge[k]) in bond_Adjacency[TrajIn.topology.atom(MultipEdge[j])]):
            zMatrix["A"][AtomIx] = k+1
            break
    for l in range(k-1, -1, -1):
        if (TrajIn.topology.atom(MultipEdge[l]) in bond_Adjacency[TrajIn.topology.atom(MultipEdge[k])]):
            zMatrix["T"][AtomIx] = l+1
            break
MultipleAtomIx = AtomIx

##  Iterate over SingleEdge Atoms   ##
for AtomIx in range(len(SingleEdge)):
    currAtom = AtomIx + MultipleAtomIx + 1

    if (TrajIn.topology.atom(MultipEdge[0]) in list(bond_Adjacency[TrajIn.topology.atom(SingleEdge[currAtom])])):
        zMatrix.loc[currAtom] = [currAtom+1, 1, 2, 3]
    elif (TrajIn.topology.atom(MultipEdge[1]) in list(bond_Adjacency[TrajIn.topology.atom(SingleEdge[currAtom])])):
        zMatrix.loc[currAtom] = [currAtom + 1, 2, 3, 4]
    else:
        #neighbor = MultipEdge_val.index(list(bond_Adjacency[TrajIn.topology.atom(SingleEdge[currAtom])])[0].index)
        print (list(bond_Adjacency[TrajIn.topology.atom(SingleEdge[currAtom])])[0])
        neighbor = MultipEdge_val.index(list(bond_Adjacency[TrajIn.topology.atom(SingleEdge[currAtom])])[0].index)
        print(neighbor+1)
        zMatrix["i"][currAtom] = currAtom+1
        zMatrix["B"][currAtom] = neighbor+1
        zMatrix["A"][currAtom] = zMatrix["B"][neighbor]
        zMatrix["T"][currAtom] = zMatrix["A"][neighbor]

print (zMatrix)

##  Turn from our indices to MDTraj atom indices
print(MultipEdge_val)
print (MultipEdge)
print (SingleEdge)
for i in ["i", "B", "A", "T"]:
    for j in range(TrajIn.n_atoms):
        if (zMatrix[i][j] is not None):
            zMatrix[i][j] = AllEdges[zMatrix[i][j]-1]

print (zMatrix)
## Generate a list of the above DoFs
DofIndices = []
##  Add Bonds   ##
for i in range(zMatrix.shape[0]):
    if (zMatrix["B"][i] is not None):
        DofIndices.append([zMatrix["i"][i],zMatrix["B"][i]])
    if (zMatrix["A"][i] is not None):
        DofIndices.append([zMatrix["i"][i],zMatrix["B"][i], zMatrix["A"][i]])
    if (zMatrix["T"][i] is not None):
        DofIndices.append([zMatrix["i"][i],zMatrix["B"][i], zMatrix["A"][i], zMatrix["T"][i]])

print (DofIndices)


## Generate an array, of shape (N_Frames,3N-6)  ##
#DoFs = np.array(TrajIn.n_frames, 3*(TrajIn.n_atoms)-6)
