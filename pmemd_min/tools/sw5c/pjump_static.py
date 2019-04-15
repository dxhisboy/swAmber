import re
import hashlib
import sys
#        bsr     $26,puts                        # [4] puts
BSR_RE = re.compile("\\s*bsr\\s+\\$26,(?P<symbol>\\S+)\\s*.*")
SYMBOL_RE = re.compile("(?P<symbol>\\S+):.*")
def process_pre_assembly(path_in, path_out = None, verbose=False):
    try:
        if verbose:
            print('pjump_static: Postprocessing %s' % path_in)
        if path_out is None:
            path_out = path_in
        
        f = open(path_in, "r")
        lines = f.readlines()
        f.close()
        f = open(path_out, "w")
        local_symbols = set([])
        for line in lines:
            lm = SYMBOL_RE.match(line)
            if lm:
                local_symbols.add(lm.groupdict()['symbol'])            
        for lineno, line in enumerate(lines):
            rm = BSR_RE.match(line)
            if rm:
                symbol = rm.groupdict()['symbol']
                if symbol not in local_symbols:
                    refhash = hashlib.md5(("%s:%d" % (path_in, lineno)).encode()).hexdigest()
                    if verbose:
                        print('pjump_static: Found a bsr reference to %s' % symbol, file=sys.stderr)
                    print("staticref.%s.hashref.%s:" % (symbol, refhash), file=f)
                    print("\tbsr $26,%s" % symbol, file=f)
                    print("\tnop", file=f)
                    print("\tnop", file=f)
                    print("\tnop", file=f)
                    print("\tnop", file=f)
                    print("\tcall $26,($27),%s" % symbol, file=f)
                else:
                    print(line.rstrip(), file=f)
            else:
                print(line.rstrip(), file=f)
        f.close()
        return True
    except:
        return False
import elf
import struct
STATICREF_RE = re.compile('staticref\\.(?P<symbol>.*)\\.hashref\\.[0-9a-f]*')
SLL27 = bytearray(struct.pack('>I', 0x1b09644b))
LDI = 0x3e
LDIH = 0x3f
def gen_mem(opc, ra, rb, disp):
    ret = opc << 26 | ra << 21 | rb << 16 | disp & 0xffff
    return bytearray(struct.pack('<I', ret))
def depart_long(x):
    parts = list(struct.unpack('HHHH', struct.pack('Q', x)))
    for i in range(3):
        if parts[i] > 0x7fff:
            parts[i] -= 0x10000
            parts[i + 1] += 1
    parts[3] &= 0xffff
    return parts
    #x0 = x & 0xffff
def process_post_assembly(path_in, path_out, verbose=True):
    try:
        print(path_in)
        elf_in = elf.ELF64(open(path_in, "rb").read())
        #print(elf_in.global_syms)
        text1_offset = elf_in.shdrs['.text1'].sh_offset
        text1_addr = elf_in.shdrs['.text1'].sh_addr
        for sym in elf_in.syms:
            if STATICREF_RE.match(sym.st_name):
                symbol = "slave_" + STATICREF_RE.match(sym.st_name).groupdict()['symbol']
                #print("%x: %s" % (sym.st_value, symbol))
                #print(symbol - test1_addr
                inst_off = (sym.st_value - text1_addr + text1_offset)
        
                if symbol in elf_in.global_syms:
                    print("%s: %s is resolved as %x" % (sym.st_name, symbol, elf_in.global_syms[symbol].st_value), file=sys.stderr)
                    parts = depart_long(elf_in.global_syms[symbol].st_value)
                    elf_in.bin[inst_off +  0 : inst_off +  4] = gen_mem(LDI , 27, 31, parts[2])
                    elf_in.bin[inst_off +  4 : inst_off +  8] = gen_mem(LDIH, 27, 27, parts[3])
                    elf_in.bin[inst_off +  8 : inst_off + 12] = SLL27
                    elf_in.bin[inst_off + 12 : inst_off + 16] = gen_mem(LDI , 27, 27, parts[0])
                    elf_in.bin[inst_off + 16 : inst_off + 20] = gen_mem(LDIH, 27, 27, parts[1])
                else:
                    print("unresolved symbol: %s" % symbol, file=sys.stderr)
                    return False
        open(path_out, "wb").write(elf_in.bin)
        return True
    except:
        return False
if __name__ == "__main__":
    #process_file("slave.s")
    process_post_assembly("./a.out", "./b.out")
