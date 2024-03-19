from NW_Part1 import Alignment

"""
Running this script will compute a pairwise alignment between two given 
sequences X and Y. 

As you hopefully see when running this script, the score of an alignment is 
returned for each pair of inputs using the method score_alignment

After confirming that running this code indeed provides a visualisation of 
a sequence alignment along with a score, begin editing this script to return an 
array of score relating to each sequence found within one of the files in ./sequences

"""

# Here is a sample output of working with the NW algorithm from Part1
x = "GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKL"
y = "NNPELQAHAGKVFKLVYEAAIQLQVTGVVVTDATLKNLGSVHVSKG"

# Here we are calling the alignment method
a = Alignment(x, y) # instantiates the alignment object, providing sequences X and Y
a.compute_alignment() # run to do the actual computation
print("Computed alignment:")
a.display_alignment() # this can be displayed to help visualise what is going on

# a.score_alignment provides a float score for X and Y alignment
print("Score of alignment: " + str(a.score_alignment()))

