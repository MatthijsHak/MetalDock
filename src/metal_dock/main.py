#!/usr/bin/env python3 

import argparse

from parser_metal_dock import Parser
from docking import docking

from monte_carlo import optimize_MC

def main():
    parser = argparse.ArgumentParser(description="Docking of organometallic compounds")
    parser.add_argument("-i","--input_file", help="Input file of .ini format specifying all parameters", required=True)

    args = parser.parse_args()

    par = Parser(args.input_file)

    if par.method.lower() == 'dock':
        docking(par)
        
    elif par.method.lower() == 'mc':
        optimize_MC(par)

    else:
        print("SPECIFY ONE OF THE TWO OPTIONS FOR MetalDock")
        print("(1) dock")
        print("(2) MC")

if __name__== '__main__':
    main()
