   .set noat
   .set noreorder
   .arch sw3
   .section .rodata, "a", "progbits"
   .align 4
   .section .text1, "ax", "progbits"
   .align 4
   .section .rodata
   .align 4
   .section .tdata_local_fix, "waT", "progbits"
   .align 0
   .section .text1, "ax", "progbits"
   .align 4
   .section .tdata_local_fix
   .org 0x0
   .align 0
   .globl constant_table_fast
   .size constant_table_fast,128
   .align 0
   .type constant_table_fast, @object
constant_table_fast:
   .quad 0xBFF71547652B82FE                           # -1.4426950408889634
   .quad 0x3FF0000000000000                           # 1.0
   .quad 0x3FDE1C2E33EFF195                           # 0.47047

   .quad 0x3FE7EE6EDE042CB5                         # 0.7478556 
   .quad 0xBFB88B941D58D65B                           # -0.0958798
   .quad 0x3FD646074B4DDC2C                           # 0.3480242

   .quad 0x3FE62E42FEFA39EF                           # 0.6931471805599453

   .quad 0x3F9E119F64B4B57F                           # 0.0293641
   .quad 0x3FDFE406A535AF0E                           # 0.4982926
   .quad 0xBFC46B957A502836                           # -0.1595332
   .quad 0xBFEFFEEC03C9B91C                           # -0.9998684

   .quad 0xBFF20DD750429B6D                           #
   .section .text1
   .align 4
   .globl simd_fast_calc_switch
   .ent simd_fast_calc_switch

simd_fast_calc_switch:
   ldi $0 , constant_table_fast($31) !tprello
   vmuld $16 , $16 , $19                              # x_square
   ldde $1 , 0($0)                                    # ln2inv
   vcpys $31 , $16 , $20                              # absx
   vmuld $19 , $1 ,  $1                               # xln2inv  
   ldde $7 , 8($0)                                    # 1.0
   
   vextf $1 , 0x1 , $2                                
   vextf $1 , 0x2 , $3
   vextf $1 , 0x3 , $4
 
   ldde $8 , 16($0)                                   # erf_q

   fcvtdlr $1 , 7 , $1
   fcvtdlr $2 , 7 , $2
   fcvtdlr $3 , 7 , $3                             
   fcvtdlr $4 , 7 , $4

   vmad $8 , $20 , $7 , $8                            # erf_t

   vshff $2 , $1 , 0x0 , $5                          
   vshff $4 , $3 , 0x0 , $6

   fcvtld $1 , $1 
   fcvtld $2 , $2
   fcvtld $3 , $3 
   fcvtld $4 , $4
   
   vshff $2 , $1 , 0x0 , $1
   vshff $4 , $3 , 0x0 , $3
   vshff $6 , $5 , 0x88 , $5                          # nln_int
   vshff $3 , $1 , 0x88 , $1                          # nln_float
   

   vdivs $7 , $8 , $8                                 # erf_t
   
   ldde $2 , 24($0)                                   # erf_a3 
   ldde $3 , 32($0)                                   # erf_a2
   ldde $4 , 40($0)                                   # erf_a1

   ldde $6 , 48($0)                                   # ln2

   ldde $21 , 56($0)                                  # exp_a4
   ldde $22 , 64($0)                                  # exp_a2
   ldde $23 , 72($0)                                  # exp_a3
   ldde $24 , 80($0)                                  # exp_a1
  

   vmad $6 , $1 , $19 , $1                            # exp_b

   vmad $2 , $8 , $3 , $2                             # erf_poly = erf_a3*erf_t + erf_a2
   
   srlow $7 , 0x34 , $27                              
   
   vmuld $1 , $1 , $6                                  # exp_b_square

   vmad $2 , $8 , $4 , $2                             # erf_poly = erf_poly*erf_t + erf_a1 

   vaddl $27 , $5 , $5                              

   vmad $21 , $6 , $22 , $21                          # exp_poly1 = exp_a4*exp_b_square + exp_a2 
   vmad $23 , $6 , $24 , $23                          # exp_poly2 = exp_a3*exp_b_square + exp_a1

   vmuld $2 , $8 , $2                                 # erf_poly = erf_poly*erf_t

   vmad $21 , $6 , $7 , $21                           # exp_poly1 = exp_poly1*exp_b_square + 1.0
   
   sllow $5 , 0x34 , $5                               # exp_p2n

   vmad $23 , $1 , $21 , $21                          # exp_poly

   ldde $25 , 88($0)                                  # neg_twice_rsqrt_pi

   vmuld $5 , $21 , $5                                # expnxsq 

   vnmad $5 , $2 , $7 , $2                            # erfabs

  
   vmuld $25 , $5 , $3                                # d_sw_dx

   vcpys $16 , $2 , $2          

   vstd $3 , 0($18) 

   vsubd $7 , $2 , $2                                 # sw

   vstd $2 , 0($17)

   ret 
   .end simd_fast_calc_switch
