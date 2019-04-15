
#commit fc4643d4ac4132668aca3fb832b3e2cf270a4530
#    auto commit for build: 190216-221816
#delete mode 100644 e.c\n
#renamed:    f.c -> g.c
import re
import os
import subprocess
COMMIT_RE = re.compile("commit (?P<rev>[0-9a-f]+)")
AUTO_COMMIT_MSG_RE = re.compile("\\s*auto commit for build: [0-9]{6}-[0-9]{6}")
DELETE_RE = re.compile("\s*delete mode\s+[0-9]+\s+(?P<filename>.*)")
#RENAME_RE = re.compile("\s*renamed:\s+(?P<src>\S+)\s*->\s*(?P<dst>\S+)")

def get_deleted(path):
    gitdir = os.path.join(path, ".git")
    try:
        gitlog = subprocess.check_output("git --git-dir %s log" % gitdir, shell=True).decode().split('\n')
        last_commit = None
        for line in gitlog:
            if COMMIT_RE.match(line):
                last_commit = COMMIT_RE.match(line).groupdict()["rev"]
            if AUTO_COMMIT_MSG_RE.match(line):
                break
        print("Last auto commit is: %s" % last_commit)
        
        gitdiff = subprocess.check_output("git --git-dir %s diff %s --summary" % (gitdir, last_commit), shell=True).decode().split('\n')
        deleted = []
        for line in gitdiff:
            if DELETE_RE.match(line):
                deleted.append(DELETE_RE.match(line).groupdict()['filename'])
            #if RENAME_RE.match(line):
            #    print("rm %s", RENAME_RE.match(line).groupdict()['src'])
        return deleted
    except:
        return []

DEPEND_RE = re.compile("\\s*(?P<obj>.*\\.o)\\s*:(?P<deps>.*)")
SPACE_RE = re.compile("\\s+")
def get_objects(path, sources):
    objects = []
    if os.path.exists(os.path.join(path, "Depends")):
        src2obj = {}
        depends = open(os.path.join(path, "Depends")).readlines()
        for line in depends:
            if DEPEND_RE.match(line):
                (obj, deps) = DEPEND_RE.match(line).groups()
                obj = obj.strip()
                deps = SPACE_RE.split(deps.strip())
                for dep in deps:
                    if dep in sources:
                        objects.append(obj)
    return objects

