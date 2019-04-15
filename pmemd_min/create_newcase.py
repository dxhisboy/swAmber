#!/usr/sw-cluster/apps/Anaconda/anaconda3/bin/python3
import sys
import os
import argparse
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools"))
import casetool

parser = argparse.ArgumentParser(description='Create a new case with specified name')
parser.add_argument("casename", action="store", help="Name of case", type=str)
parser.add_argument("-p", "--proto", action="store", help="Create case as a copy of prototype case", default=None)
args = parser.parse_args(sys.argv[1:])
#print(args)
case_vars = casetool.gen_case_vars(args.casename)
for k, v in case_vars.items():
    print("%s=%s" % (k.upper(), v))
casetool.create_case(case_vars, args.proto)

