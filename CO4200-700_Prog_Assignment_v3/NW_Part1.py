import sys
import blosum as bl

"""
Needleman-Wunsch algorithm for global alignment.
"""
class Alignment(object):
    x = None # first sequence
    y = None # second sequence
    xa = None # will be te first sequence after inserting gaps
    ya = None # and the second sequence after inserting gaps

    d = 8 # gap penalty
    columns = 80 # this no. cols for displaying the original output in shown part2

    # BLOSUM50 substitution matrix from blosum package
    blosum_50 = bl.BLOSUM(50, default=0)

    # __init__ method that sets up the value for variables used in the class
    def __init__(self, s1, s2):
        self.x = s1 # stores s1 in X
        self.y = s2 # stores s2 in Y
        self.xa = None # set xa to null
        self.ya = None # setya to null

    # substitution score for character pairs
    def s(self, x, y):
        return Alignment.blosum_50[x][y]

    # gap penalty function
    @staticmethod
    def gamma(g):
        return -g * Alignment.d

    # checks whether Xa and Ya are non-null and have equal length - basically to
    # ensure alignment has been computed first
    def check_alignment(self):
        if self.xa is None or self.ya is None:
            print("Error: alignment not yet computed")
            exit(1)
        elif len(self.xa) != len(self.ya):
            print("Error: aligned strings have different length")
            exit(1)

    # displays the alignment on the screen
    def display_alignment(self):
        self.check_alignment() # checks if Xa and Ya define an alignment

        n = len(self.xa)
        mid = [' '] * n
        i = 0
        while i < n:
            x = self.xa[i]
            y = self.ya[i]
            if x == '-' or y == '-':
                mid[i] = '.'
            elif x == y:
                mid[i] = x
            elif self.s(x, y) >= 0:
                mid[i] = '+'
            else:
                mid[i] = '.'
            i += 1

        mid_s = ''.join(mid)
        i = 0
        # output columns chars per line
        while i < n:
            endcol = i + Alignment.columns
            if endcol > n:
                endcol = n
            output_xa = ''
            output_ya = ''
            print(output_xa.join(self.xa[i:endcol]))
            print(mid_s[i:endcol])
            print(output_ya.join(self.ya[i:endcol]))
            i += Alignment.columns

    # score alignment - calls s(x,y) for substitution scores and gamma(g) for gap penalties
    # alignment must be stored in Xa and Ya
    def score_alignment(self):
        self.check_alignment()
        n = len(self.xa)
        score = 0
        i = 0
        while i < n:
            if self.xa[i] == '-':
                g = 0
                while i < n and self.xa[i] == '-':
                    g += 1
                    i += 1
                score += Alignment.gamma(g)
            elif self.ya[i] == '-':
                g = 0
                while i < n and self.ya[i] == '-':
                    g += 1
                    i += 1
                score += Alignment.gamma(g)
            else:
                score += self.s(self.xa[i], self.ya[i])
                i += 1
        return score

    # compute alignment between sequences X and Y
    def compute_alignment(self):
        # constants for remembering T/L/D in P matrix
        t = 1
        l = 2
        d = 3

        n = len(self.x) # length of X / number of chars in X
        m = len(self.y) # length of Y / number of chars in Y

        # instantiates the matrices F and P
        f = [[0] * (m + 1) for _ in range(n + 1)]
        p = [[0] * (m + 1) for _ in range(n + 1)]

        i = 0
        while i <= n:
            j = 0
            while j <= m:
                # now compute F[i][j] as in standard Needleman-Wunsch
                if i == 0 and j == 0:
                    f[i][j] = 0
                elif i == 0:
                    f[i][j] = - j * Alignment.d
                    p[i][j] = l
                elif j == 0:
                    f[i][j] = -i * Alignment.d
                    p[i][j] = t
                else:

                    diag = f[i - 1][j - 1] + self.s(self.x[i - 1], self.y[j - 1])
                    left = f[i][j - 1] - Alignment.d
                    top = f[i - 1][j] - Alignment.d
                    if diag >= left and diag >= top:
                        f[i][j] = diag
                        p[i][j] = d
                    elif left >= diag and left >= top:
                        f[i][j] = left
                        p[i][j] = l
                    else:
                        f[i][j] = top
                        p[i][j] = t
                j += 1
            i += 1

        # traceback: alignments in xb and yb
        xb = []
        yb = []
        i = n
        j = m
        while i + j > 0:
            if p[i][j] == d:
                xb.insert(0, self.x[i - 1])
                yb.insert(0, self.y[j - 1])
                i -= 1
                j -= 1
            elif p[i][j] == t:
                xb.insert(0, self.x[i - 1])
                yb.insert(0, '-')
                i -= 1
            else:
                xb.insert(0, '-')
                yb.insert(0, self.y[j - 1])
                j -= 1
        self.xa = xb
        self.ya = yb


if __name__=='__main__':
    print('This is the Needleman-Wunsch algorithm for pairwise alignment \n'
          'To run, import this program and instantiate the alignment class, \n'
          'passing two sequence to compute_alignment, this can be displayed too. '
          ''
          '\nSomething like this... ')
    print()

    # Here is a sample output of working with the NW algorithm from Part1
    x = "GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKL"
    y = "NNPELQAHAGKVFKLVYEAAIQLQVTGVVVTDATLKNLGSVHVSKG"

    # Here we are calling the alignment method
    a = Alignment(x, y)
    a.compute_alignment()
    print("Computed alignment:")
    a.display_alignment()
    print("Score of alignment: " + str(a.score_alignment()))

    # if this runs, great!
    # you will continue to part 2 of the assignment in ./Part2.py
