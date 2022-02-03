import mdtraj as md
import numpy as np
import pandas as pd

## LOAD TRAJECTORY ##
OutFN = ["./ala1_DoFs.csv", "./ala1_covmat.csv"]
TrajIn = md.load("ala1/ala1_Prod_Restraint_Global.dcd", top="ala1/ala1.prmtop")
print (TrajIn)
bondgraph = TrajIn.topology.to_bondgraph()
bond_Adjacency = bondgraph.adj

## BONDS INDICES GET ##
bonds_redundant = []
for node in bond_Adjacency:
    for firstNeighbor in bond_Adjacency[node]:
        bonds_redundant.append([node.index, firstNeighbor.index])
bonds_nonRedundant = []
for ix in range(len(bonds_redundant)):
    repeat = 0
    for jx in range(ix + 1, len(bonds_redundant)):
        if (sorted(bonds_redundant[ix]) == sorted(bonds_redundant[jx])):
            repeat = 1
            break
    if repeat == 0:
        bonds_nonRedundant.append(bonds_redundant[ix])
## BONDS INDICES END ##

print ("{} bonds found!".format(len(bonds_nonRedundant)))

## BOND LENGTHS GET ##
bonds_distances = md.compute_distances(TrajIn, bonds_nonRedundant)
## BOND LENGTHS END ##

## ANGLES INDICES GET ##
angles_redundant = []
for node in bond_Adjacency:
    for firstNeighbor in bond_Adjacency[node]:
        if (len(bond_Adjacency[firstNeighbor]) > 1):
            for secondNeighbor in bond_Adjacency[firstNeighbor]:
                if (secondNeighbor != node):
                    angles_redundant.append([node.index, firstNeighbor.index, secondNeighbor.index])
angles_nonRedundant = []
for ix in range(len(angles_redundant)):
    repeat = 0
    for jx in range(ix + 1, len(angles_redundant)):
        if sorted(angles_redundant[ix]) == sorted(angles_redundant[jx]):
            repeat = 1
            break
    if repeat == 0:
        angles_nonRedundant.append(angles_redundant[ix])
## ANGLES INDICES END ##

print ("{} angles found!".format(len(angles_nonRedundant)))

## ANGLE VALUES GET ##
angles_values = md.compute_angles(TrajIn, angles_nonRedundant)
## ANGLE VALUES END ##


## DIHEDRALS INDICES GET ##
dihedrals_redundant = []
for node in bond_Adjacency:
    for firstNeighbor in bond_Adjacency[node]:
        if (len(bond_Adjacency[firstNeighbor]) > 1):
            for secondNeighbor in bond_Adjacency[firstNeighbor]:
                if (secondNeighbor != node):
                    if (len(bond_Adjacency[secondNeighbor]) > 1):
                        for thirdNeighbor in bond_Adjacency[secondNeighbor]:
                            if (thirdNeighbor != firstNeighbor):
                                dihedrals_redundant.append([node.index, firstNeighbor.index,
                                                            secondNeighbor.index, thirdNeighbor.index])
                                # print([node.index, firstNeighbor.index, secondNeighbor.index, thirdNeighbor.index])

dihedrals_nonRedundant = []
for ix in range(len(dihedrals_redundant)):
    repeat = 0
    for jx in range(ix + 1, len(dihedrals_redundant)):
        if (sorted(dihedrals_redundant[ix]) == sorted(dihedrals_redundant[jx])):
            repeat = 1
            break
    if repeat == 0:
        dihedrals_nonRedundant.append(dihedrals_redundant[ix])
## DIHEDRALS INDICES END ##

## DIHEDRALS VALUES GET ##
dihedrals_values = md.compute_dihedrals(TrajIn, dihedrals_nonRedundant)
## DIHEDRALS VALUES END ##

print ("{} dihedrals found!".format(len(dihedrals_nonRedundant)))

## ALL DoFs CONCATENATED GET ##
All_DoFs = np.hstack((np.hstack((bonds_distances, angles_values)), dihedrals_values)).T
All_DoFs_Correlation = np.corrcoef(All_DoFs)
## Generate string indices for bonds/angles/dihedrals
for i in range(len(bonds_nonRedundant)):
    bonds_nonRedundant[i] = str(bonds_nonRedundant[i])
for i in range(len(angles_nonRedundant)):
    angles_nonRedundant[i] = str(angles_nonRedundant[i])
for i in range(len(dihedrals_nonRedundant)):
    dihedrals_nonRedundant[i] = str(dihedrals_nonRedundant[i])

DoFsIndices = bonds_nonRedundant+angles_nonRedundant+dihedrals_nonRedundant

All_DoFs_Correlation = pd.DataFrame(All_DoFs_Correlation,
                                    index = DoFsIndices,
                                    columns = DoFsIndices)
## ALL DoFs CONCATENATED END ##

## SAVE FIRST FRAME AS BAT COORDINATES ##
All_DoFs = pd.DataFrame(All_DoFs[:,0],
                        index = DoFsIndices,
                        columns = ["Value"])

## SAVE DoF Values to CSV
All_DoFs.to_csv(OutFN[0], sep="\t")
print ("DoF file (for frame 0) saved: {}".format(OutFN[0]))

## Since we're looking for corellated motion, we discard the sign
## of the covariance; i.e. a covariance of -1 is taken into account
## as much as a covariance of 1
All_DoFs_Correlation = np.abs(All_DoFs_Correlation)

## Normalize the Covariance Distribution
# CovarMin = np.min(All_DoFs_Correlation)
# CovarMax = np.max(All_DoFs_Correlation)
# All_DoFs_Correlation = (All_DoFs_Correlation - CovarMin)/(CovarMax - CovarMin)

## Return values belogning to 90th percentile
#print (All_DoFs_Correlation[All_DoFs_Correlation.index[0]])
#print (np.percentile(All_DoFs_Correlation, 99))
#print (All_DoFs_Correlation.columns[0])
for i in range(len(All_DoFs_Correlation.columns)):
    for j in range(i+1, len(All_DoFs_Correlation.index)):
        #if (All_DoFs_Correlation[All_DoFs_Correlation.columns[i]][All_DoFs_Correlation.index[j]] >= np.percentile(All_DoFs_Correlation, 0.99)):
        if (All_DoFs_Correlation[All_DoFs_Correlation.columns[i]][All_DoFs_Correlation.index[j]] >= 0.5) :
            print (All_DoFs_Correlation[All_DoFs_Correlation.columns[i]][All_DoFs_Correlation.index[j]])
            print (All_DoFs_Correlation.columns[i], All_DoFs_Correlation.index[j])

print()
## SAVE Covariance Matrix to CSV
#All_DoFs_Correlation.to_csv(OutFN[1])