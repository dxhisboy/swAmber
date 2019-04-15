#!/usr/sw-cluster/apps/Anaconda/anaconda3/bin/python3
import os
import case
import sys
import argparse
import json
import datetime
import string
sys.path.append(case.CASETOOL)
import f90nml
import runargs
import tee
SEP = os.path.sep
EPILOGS = [
    "Available datasets:",
]
datasets = os.listdir(case.INPUTROOT)
for dataset in datasets:
    EPILOGS.append("    %s" % dataset)
EPILOG = "\n".join(EPILOGS)
RUNID=datetime.datetime.strftime(datetime.datetime.now(), "%y%m%d-%H%M%S")

parser = argparse.ArgumentParser(usage="Execute the pmemd you've built", epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--dataset", help="Input dataset to use.", action="store", default="solvated_HMR")
parser.add_argument("-n", "--nprocs", help="Number of processes to use.", action="store", default="64")
parser.add_argument("-r", "--runname", help="Name of this run.", action="store", default=None)
parser.add_argument("-v", "--nlvar", metavar="'nlpath=val'", help="Temporary namelist variable (e.g. 'cntrl.nstlim=100')", action="append")
parser.add_argument("-I", "--interactive", help="Add this option for interactive run.", action="store_true")
parser.add_argument("-s", "--swlu", help="Enable SWLU.", action="store_true")
parser.add_argument("--bsubarg", help="Additional bsub arguments", action="store", type=str, default=None)
# parser.add_argument("-R", "--randomize", help="Do not use constant seeds.", action="append_const", dest="nlvar",
#                     const="cntrl.ig=-1")
parser.add_argument("-t", "--timesteps", help="Number of timesteps.", action="append", dest="nlvar",
                    type=lambda x: "cntrl.nstlim=" + str(x))

parser.set_defaults(nlvar = [])
args = parser.parse_args(sys.argv[1:])
if args.runname is None:
    args.runname = args.dataset

rundir = os.path.join(case.CASEROOT, 'run', args.runname, RUNID)
runlog = os.path.join(rundir, 'submission.log')
if not os.path.exists(rundir):
    os.makedirs(rundir)
tee.tee(runlog, sys.stdout, sys.stderr)

lastrundir = os.path.join(case.CASEROOT, 'run', args.runname, 'latest')
dataset_matched = []
for dataset in datasets:
    if dataset.startswith(args.dataset):
        dataset_matched.append(dataset)
if len(dataset_matched) < 1:
    print("No dataset matched for %s" % args.dataset)
    sys.exit(1)
if len(dataset_matched) > 1:
    print("Multiple dataset matched for %s, available datasets matched: %s" % (args.dataset, str(dataset_matched)))
    sys.exit(1)
indir = os.path.join(case.INPUTROOT, dataset_matched[0])

if os.path.lexists(lastrundir):
    os.remove(lastrundir)
os.symlink(rundir, lastrundir)
print("Run ID is: %s" % RUNID)
print("Updating namelist...")
nl = f90nml.read(indir + SEP + 'md.in')

for var in args.nlvar:
    sp = var.split('=', 1)
    if len(sp) != 2:
        print("Not recognized nlvar: " + var)
        sys.exit(1)
    path = list(map(str.strip, sp[0].split(".")))
    if len(path) < 1:
        print("Not recognized nlvar: " + var)
        sys.exit(1)        
    node = nl
    for part in path[:-1]:
        if part not in node:
            node[part] = f90nml.Namelist()
        node = node[part]
    node[path[-1]] = eval(sp[1])
os.chdir(rundir)
nl.write("md.in")

pme_args = ["-O",
            "-i", "md.in",
            "-o", "mdout",
            "-inf", "mdinfo",
            "-r", "md.rst",
            "-x", "md.nc",
            "-l", "logfile",
            "-e", "mden"]
if args.swlu:
    pme_args.append("-swlu")
dataset_json = json.JSONDecoder().decode(open(indir + SEP + "dataset.json").read())
dataset_vars = {
    "dataset" : indir
    }
for k, v in dataset_json.items():
    if k in runargs.name2arg:
        pme_args.append(runargs.name2arg[k])
        pme_args.append(string.Template(v).substitute(dataset_vars))
#print(pme_args)
#bsub -o run.log -I -J xduan_case -q q_sw_share -debug  -host_stack 256 -share_size 6144 -priv_size 4 -n 64 -b -cgsp 64  ../bin/pmemd.MPI   -O   -i md.in -o mdout_test -p solvated_HMR.prmtop -c test.rst -r md.rst -x md.nc
bsub_args = ["-debug",
             "-o", "run.log",
             "-J", os.path.basename(case.CASEROOT),
             "-q", case.QUEUE,
             "-host_stack", "256",
             "-share_size", "6144",
             "-priv_size", "4",
             "-n", str(args.nprocs),
             "-b",
             "-cgsp", "64"]
if args.bsubarg:
    bsub_args.append(args.bsubarg)
if args.interactive:
    bsub_args.append("-I")
#print(bsub_args)
pmemd_bin = case.CASEROOT + SEP + "bld" + SEP + "pmemd.MPI"
cmd = "bsub " + " ".join(bsub_args) + " " + pmemd_bin + " " + " ".join(pme_args)
print("Submitting job with: %s" % cmd)
os.system(cmd)
