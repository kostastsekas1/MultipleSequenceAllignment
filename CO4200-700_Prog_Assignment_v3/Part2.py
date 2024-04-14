import numpy as np
from NW_Part1 import Alignment


#### 2
# reading  and initializing a list to hold the sequences and adding them in a list
ListOfSequences = []
multipleSequencesFile = open("./sequences/multiple3.txt", "r")
#first line of the file is how many sequences are there
sequencesNumber = multipleSequencesFile.readline().strip()

if sequencesNumber.isnumeric():
    sequencesNumber = int(sequencesNumber)

# read the rest of file  and append sequences in the list
for i in range(0, sequencesNumber+1):
    x = multipleSequencesFile.readline().strip()
    if x != "":
        ListOfSequences.append(x)

#close file to save resources
multipleSequencesFile.close()

# function returns a matrix that when taking a list of sequences returns
# a matrix of their scores.sequence i and j can be found in the Matrix in position (i,j)
# computes the alignment for any pair of sequences even with itself
def GetScoresForPairs(SequenceList):

    tempMatrix = []
    for sequence1 in SequenceList:
        TempList = []
        for sequence2 in SequenceList:

            pair = Alignment(sequence1, sequence2)
            pair.compute_alignment()
            TempList.append(pair.score_alignment())
        tempMatrix.append(TempList)
    return tempMatrix

#calling the function and saving it in variable
ScoreMatrix = GetScoresForPairs(ListOfSequences)
# making the list an array to print it with its build in function
print(np.array(ScoreMatrix))








