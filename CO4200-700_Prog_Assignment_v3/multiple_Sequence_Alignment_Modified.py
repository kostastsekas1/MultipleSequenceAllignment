import math
import sys
import blosum as bl
import numpy as np
from sklearn.cluster import AgglomerativeClustering


class MultipleAlignment(object):
    d = 8 # gap penalty factor
    blosum_50 = bl.BLOSUM(50, default=0) # imported BLOSUM50 substitution matrix that takes char inputs x and y

    # outputs error and quits running. Used later to handle incorrect running
    @staticmethod
    def error(msg) :
        print("Error: " + msg)
        exit(1)

    # substitution score for character pairs
    # e.g. will return -d if one of the two chars is -, returns 0 if both are -,
    # and returns BLOSUM50 score otherwise
    def s(self, x, y):
        if x == '-' and y == '-':
            return 0
        if x == '-' or y == '-':
            return -MultipleAlignment.d
        if x < 'A' or x > 'Z':
            MultipleAlignment.error("s(x,y) called for illegal character x")
        if y < 'A' or y > 'Z':
            MultipleAlignment.error("s(x,y) called for illegal character y")

        return MultipleAlignment.blosum_50[x][y]

    # gap penalty
    def gamma(self, g):
        return -g * MultipleAlignment.d

    # display alignment on screen
    def displayAlignment(self, A):
        n = len(A[0])
        columns = 76
        i = 0

        # outputs columns chars per line
        while i < n:
            endcol = i + columns
            if endcol > n:
                endcol = n
            j = 0
            while j < len(A):
                print(A[j][i:endcol])
                j += 1
            print()
            i += columns

    # scores a multiple alignment, using sum-of-pairs scoring with BLOSUM50 matrix
    # and linear gap penalty
    #
    # alignment need to passed as an array of equal length
    def scoreMultipleAlignment(self, A):
        n = len(A[0])
        score = 0
        i = 0
        while i < n:
            j = 0
            while j < len(A):
                k = j + 1
                while k < len(A):
                    score += self.s(A[j][i], A[k][i])
                    k += 1
                j += 1
            i += 1
        return score

    # compute alignment between profile X and Y
    #
    # X and Y are arrays of strings
    #
    # optimal alignment of profile X to profile Y is returned as an array of strings
    def computeProfileAlignment(self, X,  Y):
        if isinstance(X, list) and isinstance(Y, list):
            n = len(X[0])-1
            m = len(Y[0])-1
            # number of strings in x
            nX = len(X)
            nY = len(Y)
            F = [[0] * (m + 1) for _ in range(n + 1)]
            P = [[0] * (m + 1) for _ in range(n + 1)]
            i = 0
            while i <= n:
                j = 0
                while j <= m:
                    # it first computes scores of column i-1 of X and
                    # column j-1 of Y with each other and with gap columns
                    scoreXtogap = 0
                    scoreYtogap = 0
                    scoreXtoY = 0
                    if i > 0:
                        k = 0
                        while k < len(X):
                            scoreXtogap += self.s(X[k][i - 1], '-') * len(Y)
                            k += 1
                    if j > 0:
                        k = 0
                        while k < len(Y):
                            scoreYtogap += self.s(Y[k][j - 1], '-') * len(X)
                            k += 1
                    if i > 0 and j > 0:
                        k = 0
                        while k < len(X):
                            l = 0
                            while l < len(Y):
                                x = X[k][i - 1]
                                y = Y[l][j - 1]

                                scoreXtoY += self.s(x, y)
                                l += 1
                            k += 1
                    # Now computes F[i][j] as standard Needleman-Wunsch
                    if i == 0 and j == 0:
                        F[i][j] = 0
                    elif i == 0:
                        F[i][j] = F[i][j - 1] + scoreYtogap
                        P[i][j] = 2
                    elif j == 0:
                        F[i][j] = F[i - 1][j] + scoreXtogap
                        P[i][j] = 1
                    else:
                        F[i][j] = F[i - 1][j - 1] + scoreXtoY
                        if F[i - 1][j] + scoreXtogap > F[i][j]:
                            F[i][j] = F[i - 1][j] + scoreXtogap
                        if F[i][j - 1] + scoreYtogap > F[i][j]:
                            F[i][j] = F[i][j - 1] + scoreYtogap
                        P[i][j] = 0
                        # store all possible maxima by setting bits 0, 1, or 2
                        # 1 = gap in Y, 2 = gap in X, 4 = no gap
                        if F[i][j] == F[i - 1][j - 1] + scoreXtoY:
                            P[i][j] += 4
                        if F[i][j] == F[i - 1][j] + scoreXtogap:
                            P[i][j] += 1
                        if F[i][j] == F[i][j - 1] + scoreYtogap:
                            P[i][j] += 2
                    j += 1
                i += 1


            # traceback: create alignment in Xtemp and Ytemp
            Xtemp = [[' '] * (n + m) for _ in range(nX)]
            Ytemp = [[' '] * (n + m) for _ in range(nY)]
            k = n + m - 1
            i = n
            j = m
            while i + j > 0:
                if P[i][j] % 4 >= 2:
                    l = 0
                    while l < nX:
                        Xtemp[l][k] = '-'
                        l += 1
                    l = 0
                    while l < nY:
                        Ytemp[l][k] = Y[l][j - 1]
                        l += 1
                    j -= 1
                    k -= 1
                elif P[i][j] % 2 > 0:
                    l = 0
                    while l < nX:
                        Xtemp[l][k] = X[l][i - 1]
                        l += 1
                    l = 0
                    while l < nY:
                        Ytemp[l][k] = '-'
                        l += 1
                    i -= 1
                    k -= 1
                elif P[i][j] >= 4:
                    l = 0
                    while l < nX:
                        Xtemp[l][k] = X[l][i - 1]
                        l += 1
                    l = 0
                    while l < nY:
                        Ytemp[l][k] = Y[l][j - 1]
                        l += 1
                    i -= 1
                    j -= 1
                    k -= 1
            k += 1
            length = n + m - k

            # now the length of the alignment is know, create the final alignment in a_final
            Afinal = [[' '] * (length) for _ in range(nX + nY)]
            q = k
            while q < n + m:
                l = 0
                while l < nX:
                    Afinal[l][q - k] = Xtemp[l][q]
                    l += 1
                l = 0
                while l < nY:
                    Afinal[l + nX][q - k] = Ytemp[l][q]
                    l += 1
                q += 1
            A = [None] * (nX + nY)
            l = 0

            # convert to array of strings and return
            while l < nX + nY:
                A[l] = ''.join(Afinal[l])
                l += 1
            return A

        # for when there is simply a string to string comparison, when one is not an array of strings
        elif isinstance(X, str) and isinstance(Y, str):
            Xa = [None] * (1)
            Xa[0] = X
            Ya = [None] * (1)
            Ya[0] = Y
            return self.computeProfileAlignment(Xa, Ya)

        # for ehen there is a string to array comparison
        elif isinstance(X, str) and isinstance(Y, list):
            Xa = [None] * (1)
            Xa[0] = X
            return self.computeProfileAlignment(Xa, Y)

        # for when there is an array to string comparison
        elif isinstance(X, list) and isinstance(Y, str):
            Ya = [None] * (1)
            Ya[0] = Y
            return self.computeProfileAlignment(X, Ya)


    # compute multiple alignment for given set of sequences
    #
    # input sequences are passed as array of strings (can have different length)
    #
    # Resulting multiple alignment must be returned as array of strings
    # (all sequences of equal length, with gap characters inserted where appropriate)
    #
    # this is where lines are selected to align with the preceeding lines, you will likely change this
    # for part 3b

    ####Part 3i
    # This is  slightly modified version of the  kimura distance implementation
    # it uses the  compute profile alignment that is produced when aligning two sequences
    # the rest is the same  and in this version there is a try except catching if xalligned and/or yalligned are  empty.
    def Kimuradistance(self,SequenceList):
        tempMatrix = []
        for sequence1 in SequenceList:
            TempList = []
            for sequence2 in SequenceList:
                pair = self.computeProfileAlignment(sequence1, sequence2)
                xalligned = pair[0]
                yalligned = pair[1]

                positions_scored = 0
                exact_matches = 0
                distance = -1
                try:
                    for i in range(0, len(xalligned)):
                        if xalligned[i] != "-" and yalligned[i] != "-":
                            positions_scored += 1
                            if xalligned[i] == yalligned[i]:
                                exact_matches += 1

                    S = exact_matches / positions_scored
                    D = 1 - S
                    distance = -math.log(1 - D - 0.2 * (D ** 2))
                except:
                    print("Xalligned and/or Yalligned is empty resulting in a div by 0")

                TempList.append(distance)
            tempMatrix.append(TempList)
        return tempMatrix


    def computeMultipleAlignment(self, X) :

        # if the input only has one sequence then do no comparison
        if len(X) <= 1:
            return X

        DistanceMatrix = self.Kimuradistance(X)

        #### Part 3ii and 3a and 3b
        ## using the distance matrix  a guide tree is produced by the  agglomerative clustering
        model = AgglomerativeClustering(linkage='complete', metric='precomputed', distance_threshold=None)
        cluster = model.fit(np.array(DistanceMatrix))

        # the clustering children that is in the form of (0, [3, 7]).
        # Can be translated it into
        # align sequence 3 with 7  in step 0  and continue in the loop
        # this makes user that A is defined for the rest of the steps
        # if  step is other than that like (2,[6,9]) align each one progressively with the rest of the cluster
        #  finally the allignment
        for i, merge in enumerate(cluster.children_):
            if i == 0:
                A = self.computeProfileAlignment(X[merge[0]],X[merge[1]])
                continue
            for index in merge:
                if index < len(X):
                    A = self.computeProfileAlignment(A, X[index])
        return A


    # main method
    @staticmethod
    def main(Args):
        in_file = None
        inputfilename = "multiple32.txt"
        print("Reading input sequences from file " + inputfilename)

        # read file
        with open(f'./sequences/{inputfilename}') as f:
            lines = f.readlines()

        # get number of sequences in file
        number_of_sequences = int(lines[0])

        N = number_of_sequences
        n = int(N)
        print('number of sequence',n)
        X = [None] * (n)
        i = 0

        while i < n:
            X[i] = lines[i+1]
            if X[i] is None:
                MultipleAlignment.error("Unexpected end of file encountered")
            i += 1

        # create MultipleAlignment object
        a = MultipleAlignment()

        # function called to construct pairwise alignment for X sequence array
        A = a.computeMultipleAlignment(X)

        print("Computed alignment:")
        a.displayAlignment(A)
        print("Score of alignment: " + str(a.scoreMultipleAlignment(A)))


if __name__=="__main__":
    MultipleAlignment.main(sys.argv)