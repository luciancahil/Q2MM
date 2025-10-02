from OpenMMStuff.Everything import run_MM
import sys

#/home/roy/anaconda3/envs/openmm/bin/python /home/roy/Documents/Q2MM/openmm/run_mm.py 0.151 317000 0.151 317000 0.151 317000 0.151 317000 1.01 1.02 1.03 1.04 1.01 1.02 1.03 1.04 1.05

float_args = [float(n) for n in sys.argv[1:]]

print(len(float_args))
error = run_MM(float_args[0], float_args[1],float_args[2],float_args[3],float_args[4],float_args[5],float_args[6],float_args[7])
print(error)
