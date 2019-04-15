#ehdr_struct = 10sHHIQQQIHHHHHH
import struct
import collections
# typedef struct
# {
#   unsigned char e_ident[EI_NIDENT];     /* Magic number and other info */
#   Elf64_Half    e_type;                 /* Object file type */
#   Elf64_Half    e_machine;              /* Architecture */
#   Elf64_Word    e_version;              /* Object file version */
#   Elf64_Addr    e_entry;                /* Entry point virtual address */
#   Elf64_Off     e_phoff;                /* Program header table file offset */
#   Elf64_Off     e_shoff;                /* Section header table file offset */
#   Elf64_Word    e_flags;                /* Processor-specific flags */
#   Elf64_Half    e_ehsize;               /* ELF header size in bytes */
#   Elf64_Half    e_phentsize;            /* Program header table entry size */
#   Elf64_Half    e_phnum;                /* Program header table entry count */
#   Elf64_Half    e_shentsize;            /* Section header table entry size */
#   Elf64_Half    e_shnum;                /* Section header table entry count */
#   Elf64_Half    e_shstrndx;             /* Section header string table index */
# } Elf64_Ehdr;

ELF64_HDR_FMT = "10sHHIQQQIHHHHHH"
ELF64_Hdr = collections.namedtuple('ELF64_Hdr', ['e_ident', 'e_type', 'e_machine', 'e_version',
                                                 'e_entry', 'e_phoff', 'e_shoff', 'e_flags',
                                                 'e_ehsize', 'e_phentsize', 'e_phnum', 'e_shentsize',
                                                  'e_shnum', 'e_shstrndx'])
# typedef struct
# {
#   Elf64_Word    sh_name;                /* Section name (string tbl index) */
#   Elf64_Word    sh_type;                /* Section type */
#   Elf64_Xword   sh_flags;               /* Section flags */
#   Elf64_Addr    sh_addr;                /* Section virtual addr at execution */
#   Elf64_Off     sh_offset;              /* Section file offset */
#   Elf64_Xword   sh_size;                /* Section size in bytes */
#   Elf64_Word    sh_link;                /* Link to another section */
#   Elf64_Word    sh_info;                /* Additional section information */
#   Elf64_Xword   sh_addralign;           /* Section alignment */
#   Elf64_Xword   sh_entsize;             /* Entry size if section holds table */
# } Elf64_Shdr;

ELF64_SHDR_FMT = "IIQQQQIIQQ"
ELF64_Shdr = collections.namedtuple('ELF64_Shdr', ['sh_name', 'sh_type', 'sh_flags', 'sh_addr',
                                                   'sh_offset', 'sh_size', 'sh_link', 'sh_info',
                                                   'sh_addralign', 'sh_entsize'])
# typedef struct
# {
#   Elf64_Word	st_name;		/* Symbol name (string tbl index) */
#   unsigned char	st_info;		/* Symbol type and binding */
#   unsigned char st_other;		/* Symbol visibility */
#   Elf64_Section	st_shndx;		/* Section index */
#   Elf64_Addr	st_value;		/* Symbol value */
#   Elf64_Xword	st_size;		/* Symbol size */
# } Elf64_Sym;


#define ELF32_ST_BIND(val)		(((unsigned char) (val)) >> 4)
#define ELF32_ST_TYPE(val)		((val) & 0xf)
#define ELF32_ST_INFO(bind, type)	(((bind) << 4) + ((type) & 0xf))

#define ELF64_ST_BIND(val)		ELF32_ST_BIND (val)
#define ELF64_ST_TYPE(val)		ELF32_ST_TYPE (val)
#define ELF64_ST_INFO(bind, type)	ELF32_ST_INFO ((bind), (type))

ELF64_ST_BIND = lambda x: (x & 0xff) >> 4
ELF64_ST_TYPE = lambda x: x & 0xf
ELF64_ST_INFO = lambda bind, type: bind << 4 | type & 0xf

STB_LOCAL      = 0               # Local symbol */
STB_GLOBAL     = 1               # Global symbol */
STB_WEAK       = 2               # Weak symbol */
STB_NUM        = 3               # Number of defined types.  */
STB_LOOS       = 10              # Start of OS-specific */
STB_GNU_UNIQUE = 10              # Unique symbol.  */
STB_HIOS       = 12              # End of OS-specific */
STB_LOPROC     = 13              # Start of processor-specific */
STB_HIPROC     = 15              # End of processor-specific */

STT_NOTYPE     = 0               # Symbol type is unspecified */
STT_OBJECT     = 1               # Symbol is a data object */
STT_FUNC       = 2               # Symbol is a code object */
STT_SECTION    = 3               # Symbol associated with a section */
STT_FILE       = 4               # Symbol's name is file name */
STT_COMMON     = 5               # Symbol is a common data object */
STT_TLS        = 6               # Symbol is thread-local data object*/
STT_NUM        = 7               # Number of defined types.  */
STT_LOOS       = 10              # Start of OS-specific */
STT_GNU_IFUNC  = 10              # Symbol is indirect code object */
STT_HIOS       = 12              # End of OS-specific */
STT_LOPROC     = 13              # Start of processor-specific */
STT_HIPROC     = 15		# End of processor-specific */


#define STN_UNDEF	0		/* End of a chain.  */

ELF64_SYM_FMT = "IBBHQQ"
ELF64_Sym = collections.namedtuple('ELF64_Sym', ['st_name', 'st_info', 'st_other', 'st_shndx', 'st_value', 'st_size'])
class ELF64:
    def split_strtab(strtab_bin):
        strings = strtab_bin.decode().split('\0')
        stringdict = {}
        str_ptr = 0
        for i in range(len(strings)):
            stringdict[str_ptr] = strings[i]
            str_ptr += len(strings[i]) + 1
        return stringdict
    def get_section_bin(elf_bin, sec_hdr):
        return elf_bin[sec_hdr.sh_offset : sec_hdr.sh_offset + sec_hdr.sh_size]
    def __init__(self, bin):
        bin = bytearray(bin)
        ehdr = ELF64_Hdr(*struct.unpack(ELF64_HDR_FMT, bin[0:64]))
        self.bin = bin
        self.ehsize    = ehdr.e_ehsize
        self.phoff     = ehdr.e_phoff
        self.phnum     = ehdr.e_phnum
        self.phentsize = ehdr.e_phentsize
        self.shoff     = ehdr.e_shoff
        self.shnum     = ehdr.e_shnum
        self.shentsize = ehdr.e_shentsize
        # shdr = []
        # #shdr_bin = bin[self.ehdr.e_shoff:]

        # 
        # for i in range(self.shoff, shdr_end, self.shentsize):
        #     shdr.append(ELF64_Shdr(bin[i:i+64]))
        shdr_end = self.shoff + self.shnum * self.shentsize
        shdr = list(map(lambda x: ELF64_Shdr(*x), struct.iter_unpack(ELF64_SHDR_FMT, bin[self.shoff : shdr_end])))
        shstrtabhdr = shdr[ehdr.e_shstrndx]
        shstrdict = ELF64.split_strtab(ELF64.get_section_bin(bin, shstrtabhdr))
        self.shdrs = {}
        for hdr in shdr:
            self.shdrs[shstrdict[hdr.sh_name]] = hdr
        symtabhdr = self.shdrs['.symtab']
        symtabbin = ELF64.get_section_bin(bin, symtabhdr)
        syms = list(map(lambda x: ELF64_Sym(*x), struct.iter_unpack(ELF64_SYM_FMT, symtabbin)))
        strtabhdr = self.shdrs['.strtab']
        strdict = ELF64.split_strtab(ELF64.get_section_bin(bin, strtabhdr))
        #self.sym_name = {}
        self.syms = list(map(lambda sym: sym._replace(st_name = strdict[sym.st_name]), syms))
        self.global_syms = {}
        for sym in self.syms:
            if ELF64_ST_BIND(sym.st_info) in [STB_GLOBAL, STB_WEAK] and ELF64_ST_TYPE(sym.st_info) == STT_FUNC:
                self.global_syms[sym.st_name] = sym

        #self.syms = syms
#        print(shstrings)
#        print(shstrdict)
        
