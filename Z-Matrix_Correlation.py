import mdtraj as md
import sys
import numpy as np
import pandas as pd

## LOAD TRAJECTORY ##
OutFN = ["./ala1_DoFs.csv"]
OutFlex = "./ala1.flex"
TrajIn = md.load("ala1/ala1_Prod_Restraint_Global.dcd",
                 top="ala1/ala1.prmtop")
print(TrajIn)
bondgraph = TrajIn.topology.to_bondgraph()

##  Make dict of indices    ##
keylist = list(bondgraph.adj.keys())
bond_Adjacency = {}
for i in range(len(keylist)):
    atomList = list(bondgraph.adj[keylist[i]].keys())
    bond_Adjacency[i] = [x.index for x in atomList]

##  In order to write a proper flex file, we need   ##
##  to know the levels of the nodes, so we always   ##
##              go Parent -> Child                  ##
rootAtom = 0
levels = np.full(TrajIn.n_atoms, None)
levels[rootAtom] = 0
for node in bond_Adjacency:
    if node != rootAtom:
        for neighbor in bond_Adjacency[node]:
            if (levels[neighbor] is not None):
                levels[node] = levels[neighbor] + 1
#print("Levels: {}".format(levels))

## BONDS INDICES GET ##
bonds_redundant = []
for node in bond_Adjacency:
    for firstNeighbor in bond_Adjacency[node]:
        bonds_redundant.append([node, firstNeighbor])
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

print("{} bonds found!".format(len(bonds_nonRedundant)))

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
                    angles_redundant.append([node, firstNeighbor, secondNeighbor])
angles_nonRedundant = []
for ix in range(len(angles_redundant)):
    repeat = 0
    for jx in range(ix + 1, len(angles_redundant)):
        if angles_redundant[ix][2] > angles_redundant[ix][0]:
            temp = angles_redundant[ix][0]
            angles_redundant[ix][0] = angles_redundant[ix][2]
            angles_redundant[ix][2] = temp
        if (angles_redundant[ix] == angles_redundant[jx]):
            repeat = 1
            break
    if (repeat == 0):
        #print (levels[angles_redundant[ix][0]])
        #print (levels[angles_redundant[ix][2]])
        if (levels[angles_redundant[ix][0]] > levels[angles_redundant[ix][2]]):
            angles_nonRedundant.append(angles_redundant[ix])
## ANGLES INDICES END ##
print("{} angles found!".format(len(angles_nonRedundant)))

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
                                dihedrals_redundant.append([node, firstNeighbor,
                                                            secondNeighbor, thirdNeighbor])
                                # print([node.index, firstNeighbor.index, secondNeighbor.index, thirdNeighbor.index])

dihedrals_nonRedundant = []
for ix in range(len(dihedrals_redundant)):
    repeat = 0
    if (levels[dihedrals_redundant[ix][0]] > levels[dihedrals_redundant[ix][2]] and
        levels[dihedrals_redundant[ix][1]] > levels[dihedrals_redundant[ix][3]]):
        dihedrals_nonRedundant.append(dihedrals_redundant[ix])
## DIHEDRALS INDICES END ##

## DIHEDRALS VALUES GET ##
dihedrals_values = md.compute_dihedrals(TrajIn, dihedrals_nonRedundant)
## DIHEDRALS VALUES END ##

print("{} dihedrals found!".format(len(dihedrals_nonRedundant)))

## ALL DoFs CONCATENATED GET ##
All_DoFs = np.vstack((np.vstack((bonds_distances.T, angles_values.T)), dihedrals_values.T))
# print (All_DoFs)
# print (All_DoFs.shape)
# sys.exit()
All_DoFs_Correlation = np.corrcoef(All_DoFs)
## Generate string indices for bonds/angles/dihedrals
bonds_nonRedundant_str = [str(x) for x in bonds_nonRedundant]
angles_nonRedundant_str = [str(x) for x in angles_nonRedundant]
dihedrals_nonRedundant_str = [str(x) for x in dihedrals_nonRedundant]
# for i in range(len(bonds_nonRedundant)):
#     bonds_nonRedundant[i] = str(bonds_nonRedundant[i])
# for i in range(len(angles_nonRedundant)):
#     angles_nonRedundant[i] = str(angles_nonRedundant[i])
# for i in range(len(dihedrals_nonRedundant)):
#     dihedrals_nonRedundant[i] = str(dihedrals_nonRedundant[i])

DoFsIndices_list = bonds_nonRedundant + angles_nonRedundant + dihedrals_nonRedundant
DoFsIndices = bonds_nonRedundant_str + angles_nonRedundant_str + dihedrals_nonRedundant_str
# print (DoFsIndices)
All_DoFs_Correlation = pd.DataFrame(All_DoFs_Correlation, index=DoFsIndices, columns=DoFsIndices)
## ALL DoFs CONCATENATED END ##

## SAVE FIRST FRAME AS BAT COORDINATES ##
All_DoFs = pd.DataFrame(All_DoFs[:, 0],
                        index=DoFsIndices,
                        columns=["Value"])

## SAVE DoF Values to CSV
# All_DoFs.to_csv(OutFN[0], sep="\t")
# print ("DoF file (for frame 0) saved: {}".format(OutFN[0]))

## Since we're looking for corellated motion, we discard the sign
## of the covariance; i.e. a covariance of -1 is taken into account
## as much as a covariance of 1
All_DoFs_Correlation = np.abs(All_DoFs_Correlation)
All_DoFs_Correlation.to_csv(OutFN[0], sep=",", float_format="%5.3f")
print("DoF correlation saved: {}".format(OutFN[0]))

## Normalize the Covariance Distribution
# CovarMin = np.min(All_DoFs_Correlation)
# CovarMax = np.max(All_DoFs_Correlation)
# All_DoFs_Correlation = (All_DoFs_Correlation - CovarMin)/(CovarMax - CovarMin)

## Return values belogning to 90th percentile
# print (All_DoFs_Correlation[All_DoFs_Correlation.index[0]])
# print (np.percentile(All_DoFs_Correlation, 99))
# print (All_DoFs_Correlation.columns[0])
correlations = 0
free_dofs = []
free_dofs_nonred = []
for i in range(len(All_DoFs_Correlation.columns)):
    for j in range(i + 1, len(All_DoFs_Correlation.index)):
        # if (All_DoFs_Correlation[All_DoFs_Correlation.columns[i]][All_DoFs_Correlation.index[j]] >= np.percentile(All_DoFs_Correlation, 0.99)):
        if (All_DoFs_Correlation[All_DoFs_Correlation.columns[i]][All_DoFs_Correlation.index[j]] >= 0.5):
            free_dofs.append(DoFsIndices_list[i])
            # print (All_DoFs_Correlation[All_DoFs_Correlation.columns[i]][All_DoFs_Correlation.index[j]])
            # print (All_DoFs_Correlation.columns[i], All_DoFs_Correlation.index[j])
            correlations += 1
print()
for ix in range(len(free_dofs) - 1):
    if free_dofs[ix] not in free_dofs[ix + 1:]:
        free_dofs_nonred.append(free_dofs[ix])
free_dofs_nonred.append(free_dofs[-1])
print(free_dofs_nonred)

print("Total bonded observables: {}".format(All_DoFs.shape[0]))
print("Correlated bonded observables: {}".format(correlations))

##  Sort joints by movement type    ##
sliders = []
anglepins = []
pins = []
universals = []
bendstretch = []
cylinders = []
sphericals = []

for ix in range(len(free_dofs_nonred)):
    if (len(free_dofs_nonred[ix]) == 2):
        if (levels[free_dofs_nonred[ix][0]] < levels[free_dofs_nonred[ix][1]]):
            sliders.append([free_dofs_nonred[ix][0],free_dofs_nonred[ix][1]])
        else:
            sliders.append([free_dofs_nonred[ix][1],free_dofs_nonred[ix][0]])

    if (len(free_dofs_nonred[ix]) == 3):
        if (levels[free_dofs_nonred[ix][0]] < levels[free_dofs_nonred[ix][2]]):
            anglepins.append([free_dofs_nonred[ix][0],free_dofs_nonred[ix][1],free_dofs_nonred[ix][2]])
        else:   #either reverse, or both on 0 and 2 on same lvl
            anglepins.append([free_dofs_nonred[ix][2],free_dofs_nonred[ix][1],free_dofs_nonred[ix][0]])

    if (len(free_dofs_nonred[ix]) == 4):
        if (levels[free_dofs_nonred[ix][0]] < levels[free_dofs_nonred[ix][3]]):
            pins.append([free_dofs_nonred[ix][1],free_dofs_nonred[ix][2]])
        else:
            if ([free_dofs_nonred[ix][2], free_dofs_nonred[ix][1]] not in pins):
                pins.append([free_dofs_nonred[ix][2], free_dofs_nonred[ix][1]])

print ("Sliders: {}\nAnglePin:{}\nPin:{}\n".format(sliders,anglepins,pins))
##  Find combined movements ##
for bond in sliders:
    if bond in anglepins:
        if bond in pins:
            sphericals.append(bond)
        else:
            bendstretch.append(bond)
    elif bond in pins:
        cylinders.append(bond)

for bond in anglepins:
    if bond in pins:
        universals.append(bond)

for bond in sphericals:
    if bond in sliders:
        sliders.remove(bond)
    if bond in anglepins:
        anglepins.remove(bond)
    if bond in pins:
        pins.remove(bond)
for bond in bendstretch:
    if bond in sliders:
        sliders.remove(bond)
    if bond in anglepins:
        anglepins.remove(bond)
for bond in cylinders:
    if bond in sliders:
        sliders.remove(bond)
    if bond in pins:
        pins.remove(bond)
for bond in universals:
    if bond in anglepins:
        anglepins.remove(bond)
    if bond in pins:
        pins.remove(bond)

print ("Sliders: {}\nAnglePin:{}\nPin:{}\nBendStretch:{}\nUniversalM:{}\nCylinder:{}\n"
       "Spherical:{}\n".format(sliders,anglepins,pins,bendstretch,universals,cylinders,
                               sphericals))
flexFile = open(OutFlex, "w")

for bond in sliders:
    flexFile.write
for ix in range(len(free_dofs_nonred)):
    if (len(free_dofs_nonred[ix]) == 2):
        if (levels[free_dofs_nonred[ix][0]] < levels[free_dofs_nonred[ix][1]]):
            sliders.write("{} {} Slider\n".format(free_dofs_nonred[ix][0],
                                                   free_dofs_nonred[ix][1]))
        else:
            flexFile.write("{} {} Slider\n".format(free_dofs_nonred[ix][1],
                                                   free_dofs_nonred[ix][0]))

    if (len(free_dofs_nonred[ix]) == 3):
        if free_dofs_nonred[ix][1:] not in jointsAngle:
            jointsAngle.append(free_dofs_nonred[ix][1:])
            if (levels[free_dofs_nonred[ix][1]] < levels[free_dofs_nonred[ix][2]]):
                flexFile.write("{} {} AnglePin\n".format(free_dofs_nonred[ix][1],
                                                         free_dofs_nonred[ix][2]))
            else:   #either reverse, or both on 0 and 2 on same lvl
                flexFile.write("{} {} AnglePin\n".format(free_dofs_nonred[ix][2],
                                                         free_dofs_nonred[ix][1]))
    if (len(free_dofs_nonred[ix]) == 4):
        if free_dofs_nonred[ix][1:3] not in jointsPin:
            jointsPin.append(free_dofs_nonred[ix][1:3])
            if (levels[free_dofs_nonred[ix][1]] < levels[free_dofs_nonred[ix][2]]):
                flexFile.write("{} {} Pin\n".format(free_dofs_nonred[ix][1],
                                                    free_dofs_nonred[ix][2]))
            else:
                flexFile.write("{} {} Pin\n".format(free_dofs_nonred[ix][2],
                                                    free_dofs_nonred[ix][1]))
flexFile.close()
