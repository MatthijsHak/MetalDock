import argparse

from docking import docking
from train_GA import train_GA
# from test_GA import test_GA 

def main():
    parser = argparse.ArgumentParser(description="Force Field Optimization System")
    parser.add_argument("-i","--input_file", help="Forcebalance input file", required=True)

    option = parser.add_mutually_exclusive_group(required=True)
    option.add_argument("-dock", action='store_true', help="Initialize docking procedure with metal ligand.")
    option.add_argument("-train", action='store_true', help="Initialize training procedure with GA to otbain LJ parmaters for a specific metal.")
    option.add_argument("-test",  action='store_true', help="Initialize testing procedure for a specific metal.")

    args = parser.parse_args()

    if args.dock == True:
        docking(args.input_file)

    if args.train == True:
        train_GA(args.input_file)

    if args.test == True:
        test_GA(args.input_file)


if __name__== '__main__':
    main()
