import numpy as np
import string

# helper function to generate dummy testing data
def generate_background(size, unique):
    alphabet = string.ascii_uppercase
    rng = np.random.default_rng() # random generator
    indexes = rng.integers(low=0, high=len(alphabet), size=unique) # choose unique letters
    print(alphabet[indexes])

generate_background(1,2)