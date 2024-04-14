import math

import numpy as np

from NW_Part1 import Alignment

# this is  the same as for part 2
ListOfSequences = []
multipleSequencesFile = open("./sequences/multiple10.txt", "r")
sequencesNumber = multipleSequencesFile.readline().strip()

if sequencesNumber.isnumeric():
    sequencesNumber = int(sequencesNumber)

for i in range(0, sequencesNumber+1):
    x = multipleSequencesFile.readline().strip()
    if x != "":
        ListOfSequences.append(x)

multipleSequencesFile.close()

#### 3a
# Calculating the kimura distance of an allignment and creates a matrix
#  it takes in to consideration the positions_scored if there are no gaps
# and the exact_matches if a residue is the same with another
# it produces a matrix

def Kimuradistance(SequenceList):
    tempMatrix = []
    for sequence1 in SequenceList:
        TempList = []
        for sequence2 in SequenceList:
            pair = Alignment(sequence1, sequence2)
            #same as part 2 comput the alignment between pairs
            pair.compute_alignment()
            #xalligned and yalligned are the  sequences with gaps if needed from the alignment

            xalligned = pair.xa
            yalligned = pair.ya
            # position scored and exact matches are  counters  and are explained above
            positions_scored = 0
            exact_matches = 0
            for i in range(0, len(xalligned)):
                if xalligned[i] != "-" and yalligned[i] != "-":
                    positions_scored += 1
                    if xalligned[i] == yalligned[i]:
                        exact_matches += 1

            # then the rest is calculated and this is done for all the  pairs of sequences
            S = exact_matches / positions_scored
            D = 1 - S
            distance = -math.log(1-D-0.2*(D**2))

            TempList.append(distance)
        tempMatrix.append(TempList)
    return tempMatrix

# matrix is saved and printed like  part2
DistanceMatrix = Kimuradistance(ListOfSequences)
print(np.array(DistanceMatrix))
