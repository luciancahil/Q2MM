import numpy as np
import argparse
from OpenMMStuff.Everything import run_MM
import os

def square_sin(param, x):
    assert(len(x) == 2)
    return (x[0] - 4.8)**2 + 2 * (np.sin(x[1] + 0.1))**2


def open_mm_forces(x):
    assert(len(x) == 17)


    bonds = tuple(x[-9:-5])
    angles = tuple(x[-5:])
    print(len(x))
    error = run_MM(x[0], x[1],x[2],x[3],x[4],x[5],x[6],x[7], bonds, angles)
    return (error)


def open_mm_forceless(x):
    assert(len(x) == 17)


    bonds = tuple(x[-9:-5])
    angles = tuple(x[-5:])
    print(len(x))
    error = run_MM(x[0], x[1],x[2],x[3],x[4],x[5],x[6],x[7], bonds, angles, False)
    return (error)

def open_mm_focused(x):
    print(len(x))
    assert(len(x) == 8)


    bonds = tuple(x[-9:-5])
    angles = tuple(x[-5:])
    print(len(x))
    error = run_MM(x[0], x[1],x[2],x[3],x[4],x[5],x[6],x[7], use_forces=False)
    return (error)

def choose_function(name):
    if name == "square_sin":
        return square_sin
    else:
        raise(ValueError("Function not defined"))


def get_input_vals(name):
    input_file = os.path.join("BO_data", name, "next.csv")
    try:
        input_file = open(input_file)
    except(FileNotFoundError):
        input_file = os.path.join("BO_data", name, "initial.csv")
        input_file = open(input_file)   

    line = input_file.readline()

    return [float(part) for part in line.split(",")]

def append(name, x, y):
    file = os.path.join("BO_data", name, "history.csv")
    file = open(file, mode='a')

    x.append(y)

    val_list = ",".join([str(n) for n in x])
    file.write(val_list + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A simple script with arguments.")
    parser.add_argument("--name", type=str, help="name of function", required=True)
    parser.add_argument("--parameter", type=str, required=False, default=None)
    args = parser.parse_args()
    name = args.name

    if(args.parameter != "Dummy"):
        name = name + "-" + args.parameter



    x = get_input_vals(name)
    black_box_function = choose_function(args.name)


    y = black_box_function(args.parameter, x)


    print("f({}) = {}".format(x, y))
    append(name, x, y)
