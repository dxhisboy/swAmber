#include "/usr/sw-mpp/swcc/swgcc-binary-sw_64.1206/sw_64thl-unknown-linux-gnu/sys-include/asm/regdef.h"
        .set noat
        .set noreorder
        .arch sw3
        .section .rodata, "a", "progbits"
        .align 4
        .section .text1, "ax", "progbits"
        .align 4
        .section .rodata
        .align 4
        .macro do_prefetch off
        vldd t4, \off(t0)
        vldd t5, \off(t1)
        vldd t6, \off(t2)
        vldd t7, \off(t3)
        vshff t5, t4, 0x44, t8
        vshff t5, t4, 0xee, t9
        vshff t7, t6, 0x44, t10
        vshff t7, t6, 0xee, t11
        vshff t10, t8, 0x88, t4
        vshff t10, t8, 0xdd, t5
        vshff t11, t9, 0x88, t6
        vshff t11, t9, 0xdd, t7
        vstd t4, ((4 * \off) +  0)(a1)
        vstd t5, ((4 * \off) + 32)(a1)
        vstd t6, ((4 * \off) + 64)(a1)
        vstd t7, ((4 * \off) + 96)(a1)        
        .endm

        .section .text1, "ax", "progbits"
        .align 4
        .globl parameters_prefetch_32
        .ent parameters_prefetch_32

parameters_prefetch_32:
        subl a2, 1, t0
        subl a3, 1, t1
        subl a4, 1, t2
        subl a5, 1, t3
        sll t0, 8, t0 
        sll t1, 8, t1
        sll t2, 8, t2
        sll t3, 8, t3
        addl t0, a0, t0 
        addl t1, a0, t1
        addl t2, a0, t2
        addl t3, a0, t3
        do_prefetch   0
        do_prefetch  32
        do_prefetch  64
        do_prefetch  96
        do_prefetch 128
        do_prefetch 160
        do_prefetch 192
        do_prefetch 224
        ;; vinsf t0, v0, 0, v0
        ;; vinsf a2, v0, 1, v0
        ;; vinsf a0, v0, 2, v0
        ;; vinsf t0, v0, 3, v0
        ret
        .end parameters_prefetch_32
