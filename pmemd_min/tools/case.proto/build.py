#!/usr/sw-cluster/apps/Anaconda/anaconda3/bin/python3
import os
import datetime
import sys
import io
import subprocess
import case
sys.path.append(case.CASETOOL)
import filescan
import tee
#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

BLDID=datetime.datetime.strftime(datetime.datetime.now(), "%y%m%d-%H%M%S")
bldlog = os.path.join(case.CASEROOT, "bld", "bldlog.%s" % BLDID)
objdir = os.path.join(case.CASEROOT, "bld", "obj")
if not os.path.exists(objdir):
    os.makedirs(objdir)

tee.tee(bldlog, os.path.sys.stdout, sys.stderr)

os.environ["CASEROOT"] = case.CASEROOT
os.environ["PMEMDSRC"] = case.PMEMDSRC
os.environ["CASETOOL"] = case.CASETOOL
os.environ['BLDID']=BLDID
os.environ["BINDIR"] = os.path.join(case.CASEROOT, "bld")

os.chdir(case.CASEROOT)
print("Commiting SourceMods before build...")
os.chdir("SourceMods")
os.system("git add -A .")
deleted_files = filescan.get_deleted('.')

os.system("git commit --allow-empty -m 'auto commit for build: %s'" % BLDID)

os.chdir(objdir)
if len(deleted_files):
    print("File deleted in SourceMods since last build: %s" % " ".join(deleted_files))
    objects_to_delete = filescan.get_objects(objdir, deleted_files)
    print("Removing corresponding objects...")
    for objfile in objects_to_delete:
        os.remove(objfile)
        print("Removed: %s" % objfile)

filepath = open("Filepath", "w")
filepath.write(os.path.join(case.CASEROOT, "SourceMods") + "\n")
filepath.write(case.PMEMDSRC + "\n")
filepath.write(os.path.join(case.PMEMDSRC, "include") + "\n")
filepath.close()

make_cmd = "make -r -f %s/Makefile -j %d" % (case.CASETOOL, case.MAKE_J)
print("Executing: %s" % make_cmd)
make = subprocess.Popen(make_cmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
if make.wait():
    print("Build failed, EXITTING!!!")
    exit(1)
else:
    print("Done")
