import os
import sys
import string
import shutil
import config

def get_root():
    filedir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(filedir, "..")

def gen_case_vars(casename):
    root = get_root()
    casetool  = os.path.abspath(os.path.join(root, "tools"))
    caseroot  = os.path.abspath(os.path.join(root, "cases", casename))
    pmemdsrc  = os.path.abspath(os.path.join(root, "src"))
    inputroot = os.path.abspath(config.INPUTROOT)
    return {
        "pmemdroot" : root,
        "inputroot" : inputroot,
        "caseroot" : caseroot,
        "casetool" : casetool,
        "pmemdsrc" : pmemdsrc
        }
def gen_casepy(case_vars):
    caseconf_initial = '''
INPUTROOT="$inputroot"
CASEROOT="$caseroot"
CASETOOL="$casetool"
PMEMDSRC="$pmemdsrc"
MAKE_J=32
QUEUE="q_sw_share"
'''

    casepy_template = string.Template(caseconf_initial)
    casepy = casepy_template.substitute(case_vars)
    return casepy

def create_case(case_vars, proto=None):
    parent = os.path.abspath(os.path.join(case_vars["caseroot"], ".."))
    if not os.path.exists(parent):
        print("cases directory not exists, creating...")
        os.makedirs(parent)
    if os.path.exists(case_vars["caseroot"]):
        print("Case directory already exists, EXITTING!!!")
        sys.exit(1)
    if proto is None:
        print("Copying case prototype to caseroot...")
        shutil.copytree(os.path.join(case_vars["casetool"], "case.proto"), case_vars["caseroot"])
    else:
        print("Copying %s to caseroot..." % proto)
        shutil.copytree(proto, case_vars["caseroot"], ignore=shutil.ignore_patterns("run", "pmemd.MPI.*", "bldlog.*", "src", ".git"))
        print("Copying .git to SourceMods...")
        shutil.copytree(os.path.join(proto, "SourceMods", ".git"), os.path.join(case_vars["caseroot"], "SourceMods", ".git"), copy_function=shutil.copy)
    print("Generating case.py...")
    casepy = gen_casepy(case_vars)
    open(os.path.join(case_vars["caseroot"], "case.py"), "w").write(casepy)
    print("Linking src directory to case...")
    os.symlink(case_vars["pmemdsrc"], os.path.join(case_vars["caseroot"], "src"))
    print("Done")
#print(gen_case_vars("test"))
#create_case(gen_case_vars("test"))

