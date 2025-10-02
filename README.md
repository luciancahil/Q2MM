# Bayes optimization


## How to use

1. Decide on a name for your project. For example, square_sin.
1. Make 2 files [name]_history.txt and [name]_next.txt (so square_sin_history.txt)
1. In the first line of [name]_history.txt, write the dimensionality of the input (for example, 2 if the input has 2 dimensions).
1. In every subsequent line, define the range of each variable, seperated by a comma (if the range is from 0 - 6, the line should be 0,6)
1. In [name]_next.txt, give an initial input for the file. Could just be 0,0 for a 2D function.
1. In function.py, code the function you want minimized as a global function. 
1. In choose_function in function.py, add to the if statement a combination of name to return you r newly defined function
1. Run bash run_bo.sh [name] to begin running. For example, if your function name is "square_sin", type "run_bo.sh square_sin" in your terminal