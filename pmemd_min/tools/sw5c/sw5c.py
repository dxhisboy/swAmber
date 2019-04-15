#!/usr/sw-cluster/apps/Anaconda/anaconda3/bin/python3
import os
import sys
import re
import argparse
import enum
import collections
import tempfile
import pjump_static
class TargetFile(enum.Enum):
    PREPROCESSED = 1
    ASSEMBLY = 2
    OBJECT = 3
    EXECUTABLE = 4

class TargetPE(enum.Enum):
    MPE = 1
    CPE = 2
    ACC = 3
    AIO = 4

class CompilerType(enum.Enum):
    C = 1
    CPP = 2
    FORTRAN = 3
    ASM = 4
    OBJ = 5

InputFile = collections.namedtuple('InputFile', ['name', 'prefix', 'type'])
TempFile = collections.namedtuple('TempFile', ['name', 'origin'])
parser = argparse.ArgumentParser(usage='This is a wrapped compiler script.\n', allow_abbrev=False)

# target file types
parser.add_argument('-E', action='store_const', dest='target_filetype', const=TargetFile.PREPROCESSED)
parser.add_argument('-S', action='store_const', dest='target_filetype', const=TargetFile.ASSEMBLY)
parser.add_argument('-c', action='store_const', dest='target_filetype', const=TargetFile.OBJECT)
parser.set_defaults(target_filetype=TargetFile.EXECUTABLE)

# target pe types
parser.add_argument('-host',  action='store_const', dest='target_petype', const=TargetPE.MPE)
parser.add_argument('-slave', action='store_const', dest='target_petype', const=TargetPE.CPE)
parser.add_argument('-acc',   action='store_const', dest='target_petype', const=TargetPE.ACC)
parser.add_argument('-aio',   action='store_const', dest='target_petype', const=TargetPE.AIO)
parser.set_defaults(target_petype=TargetPE.MPE)

parser.add_argument('-o', action='store', dest='output')
parser.set_defaults(output=None)
# mpi disable
# parser.add_argument('-nompi', action='store_false', dest='mpi')
# parser.set_defaults(mpi=True)

# module generation path
parser.add_argument('-J', action='store', dest='module_path')
parser.set_defaults(module_path='./')

# verbose level
parser.add_argument('-v', action='store_const', const=1, dest='verbose')
parser.add_argument('-vv', action='store_const', const=2, dest='verbose')
parser.add_argument('-keep', action='store_true', dest='keep')
#parser.add_argument('-dry', action='store_true', dest='dry')
parser.set_defaults(verbose=0, dry=False)

# static pjump hook
parser.add_argument('-pjump_static', action='store_true', dest='pjump_static')
parser.set_defaults(pjump_static=False)

INPUT_RE = re.compile('(?P<prefix>.*)\\.(?P<suffix>o|c|h|cpp|cc|C|s|S|F90|f90|F|f|F77|f77)')

parsed_args, unparsed_arg = parser.parse_known_args(sys.argv[1:])
#print(parsed_args, unparsed_arg)
very_verbose = False
if parsed_args.verbose == 2:
    very_verbose = True
temp_dir = tempfile.mkdtemp(prefix='sw5mpicc')
# if not os.path.exists(temp_dir):
#     os.mkdir("mytmp")
def exit(msg = None, val = 0):
    if msg:
        print(msg, file=sys.stderr)
    if not parsed_args.keep:
        os.system("rm -r %s" % temp_dir)
    else:
        print("Temperory file leaved at %s" % temp_dir)
    merge_timestamp()
    sys.exit(val)

def execute(args):
    if parsed_args.verbose > 0 or parsed_args.dry:
        print(' '.join(args), file=sys.stderr)
    if not parsed_args.dry:
        ret = os.system(' '.join(args))
        if ret != 0:
            exit('Execution "%s" failed with exit code %d (killed by sig %d)' % (' '.join(args), ret >> 8, ret & 0xff), val=1)
output_list = []
def merge_timestamp():
    if len(output_list) > 1:
        execute(["touch -r", " ".join(output_list)])
def find_compiler(input_file):
    ret = []
    if parsed_args.target_petype == TargetPE.MPE:
        if input_file.type == CompilerType.C:
            ret.append('mpicc')
        elif input_file.type == CompilerType.CPP:
            ret.append('mpiCC')
        elif input_file.type == CompilerType.FORTRAN:
            ret.append('mpif90')
            ret.append('-J')
            ret.append(temp_dir)
    elif parsed_args.target_petype == TargetPE.CPE:
        if input_file.type == CompilerType.C:
            ret.append('sw5cc')
            ret.append('-slave')
        elif input_file.type == CompilerType.CPP:
            exit("CPE compilers do not handle c++ yet!", val = 1)
        elif input_file.type == CompilerType.FORTRAN:
            ret.append('sw5f90')
            ret.append('-slave')
    elif parsed_args.target_petype == TargetPE.ACC:
        if input_file.type == CompilerType.C:
            ret.append('swacc')
        if input_file.type == CompilerType.CPP:
            ret.append('swaCC')
        if input_file.type == CompilerType.FORTRAN:
            ret.append('swafort')
    return ret
def find_assembler():
#    ret = ['sw5cc', '-c']
    if parsed_args.target_petype == TargetPE.CPE:
#        ret.append('-slave')
#        return ['/usr/sw-mpp/swcc/sw5gcc-binary/bin/sw5as']
        return ['sw5cc', '-slave', '-c']
    else:
        return ['sw5cc', '-host', '-c']
#        return ['/usr/sw-mpp/swcc/swgcc-binary/bin/swas']
#        ret.append('-host')
    return ret

def find_assembler_asm(path):
#    ret = ['sw5cc', '-c']
    if parsed_args.target_petype == TargetPE.CPE:
#        ret.append('-slave')
        return ['/usr/sw-mpp/swcc/sw5gcc-binary/bin/sw5as']
    else:
        return ['/usr/sw-mpp/swcc/swgcc-binary/bin/swas']
#        ret.append('-host')
    return ret

# Input file find phase:
file_list = []
compiler_args = []
for unparsed_args in unparsed_arg:
    m = INPUT_RE.fullmatch(unparsed_args)
    if m:
        suffix = m.groupdict()['suffix']
        prefix = m.groupdict()['prefix']
        compiler = None
        if suffix in ['c', 'h']:
            compiler = CompilerType.C
        elif suffix in ['cpp', 'cc', 'C']:
            compiler = CompilerType.CPP
        elif suffix[0] in ['F', 'f']:
            compiler = CompilerType.FORTRAN
        elif suffix in ['s', 'S']:
            compiler = CompilerType.ASM
        elif suffix in ['o']:
            compiler = CompilerType.OBJ
        file_list.append(InputFile(unparsed_args, prefix, compiler))
    else:
        compiler_args.append(unparsed_args)
    if parsed_args.verbose == 2:
        compiler_args.append('-v')

# Preprocessing phase
if parsed_args.target_filetype == TargetFile.PREPROCESSED:
    for input_file in file_list:
        args = find_compiler(input_file)
        args.append('-E')
        args.append(input_file.name)
        args.extend(compiler_args)
        execute(args)
    exit()

# Assembly generation phase:
for input_file in file_list:
    if input_file.type in [CompilerType.C, CompilerType.CPP, CompilerType.FORTRAN]:
        #execute(['cp', input_file.name, temp_dir + os.path.sep + os.path.basename(input_file.name) + '.s'])
#        continue
        args = find_compiler(input_file)
        args.append('-S')
        args.append(input_file.name)
        args.extend(compiler_args)
        if parsed_args.target_filetype == TargetFile.ASSEMBLY:
            if parsed_args.output is not None:
                args.append('-o')
                args.append(parsed_args.output)
                output_list.append(parsed_args.output)
            else:
                args.append('-o')
                args.append(os.path.basename(input_file.name) + '.s')
                output_list.append(os.path.basename(input_file.name) + '.s')
        else:
            args.append('-o')
            args.append(temp_dir + os.path.sep + os.path.basename(input_file.name) + '.s')
        execute(args)
        if input_file.type == CompilerType.FORTRAN:
            #modules = filter(lambda path: path.endswith('.mod') and path not in modules_for_update, os.listdir(temp_dir))
            modules = filter(lambda path: path.endswith('.mod'), os.listdir(temp_dir))
            module_path = parsed_args.module_path or './'
            for mod in modules:
                mod_new = temp_dir + os.path.sep + mod
                mod_old = module_path + os.path.sep + mod
                mod_bin_new = open(mod_new, "rb").read()
                if os.path.exists(mod_old):
                    mod_bin_old = open(mod_old, "rb").read()
                    if mod_bin_new != mod_bin_old:
                        if parsed_args.verbose:
                            print("Adding module file %s for updating" % mod_old)
                        #modules_for_update[mod_old] = mod_bin_new
                        open(mod_old, "wb").write(mod_bin_new)
                        output_list.append(mod_old)
                else:
                    #modules_for_update[mod_old] = mod_bin_new
                    open(mod_old, "wb").write(mod_bin_new)
                    output_list.append(mod_old)
    elif input_file.type in [CompilerType.ASM]:
        execute(['cp', input_file.name, temp_dir + os.path.sep + os.path.basename(input_file.name)])
if parsed_args.target_filetype == TargetFile.ASSEMBLY:
    exit()

if parsed_args.pjump_static:
    for input_file in file_list:
        pjump_static.process_pre_assembly(temp_dir + os.path.sep + os.path.basename(input_file.name) + '.s', verbose=very_verbose)

# Object generation phase:
for input_file in file_list:
    if input_file.type in [CompilerType.C, CompilerType.CPP, CompilerType.FORTRAN, CompilerType.ASM]:
        args = find_assembler()
        if input_file.type == CompilerType.ASM:
            args.append(temp_dir + os.path.sep + os.path.basename(input_file.name))
        else:
            args.append(temp_dir + os.path.sep + os.path.basename(input_file.name) + '.s')
        if parsed_args.target_filetype == TargetFile.OBJECT:
            if parsed_args.output is not None:
                args.append('-o')
                args.append(parsed_args.output)
                output_list.append(parsed_args.output)
            else:
                args.append('-o')
                args.append(os.path.basename(input_file.prefix) + '.o')
                output_list.append(os.path.basename(input_file.prefix) + '.o')
        else:
            args.append('-o')
            args.append(temp_dir + os.path.sep + os.path.basename(input_file.name) + '.o')
        #args.extend(compiler_args)
        execute(args)
    elif input_file.type in [CompilerType.OBJ]:
        execute(['cp', input_file.name, temp_dir + os.path.sep + os.path.basename(input_file.name) + '.o'])
        
# # Update modules after object files generated
# for path, bin in modules_for_update.items():
#     try:
#         open(path, 'wb').write(bin)
#         output_list.append(path)
#     except:
#         exit("Failed to update module %s" % path, val=1)

if parsed_args.target_filetype == TargetFile.OBJECT:
    exit()

# Link phase:
if len(file_list) > 0:
    link_args = ['mpicc', '-hybrid']

    for input_file in file_list:
        link_args.append(temp_dir + os.path.sep + os.path.basename(input_file.name) + '.o')
    link_args.extend(compiler_args)
    if parsed_args.output:
        link_args.append('-o')
        link_args.append(parsed_args.output)
    link_args.append('-lfortran -lffio -lstdc++ -lmpifort -lfortran_slave -lffio_slave')
    execute(link_args)
    if parsed_args.pjump_static:
        pjump_static.process_post_assembly(parsed_args.output or './a.out', parsed_args.output or './a.out', verbose=very_verbose)
    output_list.append(parsed_args.output or './a.out')
else:
    #print("No input files", file=sys.stderr)
    exit("No input files", val=1)
exit(0)
