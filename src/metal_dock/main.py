#!/usr/bin/env python3 

import argparse

from parser_metal_dock import Parser
from docking import docking
from train_GA import train_GA
from test_GA import test_GA

def main():
    parser = argparse.ArgumentParser(description="Docking of organometallic compounds")
    parser.add_argument("-i","--input_file", help="Input file of .ini format specifying all parameters", required=True)

    args = parser.parse_args()

    par = Parser(args.input_file)

    if par.method.lower() == 'dock':
        docking(par)
        
    elif par.method.lower() == 'train':
        train_GA(par)

    elif par.method.lower() == 'test':
        test_GA(par)

    else:
        print("SPECIFY ONE OF THE THREE OPTIONS FOR MetalDock")
        print("(1) dock")
        print("(2) train")
        print("(3) test")

if __name__== '__main__':
    main()
