"""Script to generate the implicit (unencoded) and explicit affine layers of a CRAX instance for a fixed key."""
from collections import namedtuple
from functools import partial

import os

import sage.all

from boolcrypt.utilities import (
    substitute_variables, BooleanPolynomialRing, vector2int,
    int2vector, compose_affine, matrix2anf, compose_anf_fast, get_smart_print
)

from boolcrypt.modularaddition import get_implicit_modadd_anf

from argparse import ArgumentParser

SpeckInstance = namedtuple('SpeckInstance', 'name, default_rounds, ws, m, alpha, beta')

# CRAXInstance = namedtuple('CRAXInstance','name,default_rounds,Shrift[8],ws')

class CRAX:
    default_rounds = 50
    Shift=[31,24,17,17,0,31,24,16]
    ws = 32
    
class Pair:
    def __init__(self):
        self.l = 0
        self.r = 0
    def __init__(self, left, right):
        self.l = left
        self.r = right
   


speck_instances = {
    8: SpeckInstance("Speck_8_16", 4, 4, 4, 2, 1),  # non-standard
    32: SpeckInstance("Speck_32_64", 22, 16, 4, 7, 2),
    64: SpeckInstance("Speck_64_128", 27, 32, 4, 8, 3),
    128: SpeckInstance("Speck_128_256", 34, 64, 4, 8, 3),
}


def get_round_keys(master_key):#128bit Master_key input as a tuple of 4 elements
    default_rounds = CRAX.default_rounds
    n = CRAX.ws
    constant = [0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738, 0xBB1185EB]
    round_keys = [Pair(0,0) for _ in range(51)]
    i=0
    j=0
    while i <default_rounds:
        
        round_keys[i].l = master_key[2*(j%2)+0]^j
        round_keys[i].r = master_key[2*(j%2)+1]
        i = i+1
        round_keys[i].l = constant[j%5]
        round_keys[i].r = 0
        i = i+1
        
        round_keys[i].l = constant[j%5]
        round_keys[i].r = 0
        i = i+1
        
        round_keys[i].l = constant[j%5]
        round_keys[i].r = 0
        i = i+1
        
        round_keys[i].l = constant[j%5]
        round_keys[i].r = 0
        i = i+1
        
        j = j+1
    round_keys[i].l = master_key[0]
    round_keys[i].r = master_key[1]
    

    return round_keys


def get_implicit_unencoded_affine_layers(
        master_key, only_x_names=False,
        return_also_explicit_affine_layers=False,
        return_implicit_round_functions=False  # only needed for debugging
):
    n = CRAX.ws
    rounds = CRAX.default_rounds

    ws = n
    # generating shift operation in Matrices{Step 1:In 1 word}
    bpr = sage.all.GF(2)
    
    identity_matrix = partial(sage.all.identity_matrix, bpr)
    zero_matrix = partial(sage.all.zero_matrix, bpr)

    rightshift_values = [31,24,17,17,0,31,24,16]
    RightS_Matrices = []
    for RS_values in rightshift_values:
        matrix = sage.all.block_matrix(bpr, 2, 2, [
        [zero_matrix(ws - RS_values, RS_values), identity_matrix(ws - RS_values)],
        [identity_matrix(RS_values), zero_matrix(RS_values, ws - RS_values)]])
        RightS_Matrices.append(matrix)
    
    leftshift_values = [31,0,17,0,0,0,24,0]
    LeftS_Matrices = []
    for LS_values in leftshift_values:
        matrix = sage.all.block_matrix(bpr, 2, 2, [
        [zero_matrix(LS_values, ws - LS_values), identity_matrix(LS_values)],
        [identity_matrix(ws - LS_values), zero_matrix(ws - LS_values, LS_values)]])
        LeftS_Matrices.append(matrix)
    # generating shift operation in Matrices {Step 1...END}
    ii = identity_matrix(ws)
    zz = zero_matrix(ws, ws)

    # generating shift operation in Matrices {Step 2:In 2 word}
    # in the rounds of Alzette
    # 1st round
    identity_rotateright31_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, RightS_Matrices[0]]])#{0: 31}
    # 2nd round
    identity_rotateleft31_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, LeftS_Matrices[0]]])#{0: 31} 
    identity_rotateright24xor_Matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [RightS_Matrices[1], ii]])#{1: 24}
    identity_rotateright17_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, RightS_Matrices[2]]])#{2: 17}
    # 3rd round
    identity_rotateleft17_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, LeftS_Matrices[2]]])#{2: 17} 
    identity_rotateright17xor_Matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [RightS_Matrices[3], ii]])#{3: 17}
    identity_rotateright0_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, RightS_Matrices[4]]])#{4: 0} 
    # 4th round 
    identity_rotateleft0_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, LeftS_Matrices[4]]])#{4: 0}
    identity_rotateright31xor_Matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [RightS_Matrices[5], ii]])#{5: 31}
    identity_rotateright24_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, RightS_Matrices[6]]])#{6: 24} 
    # 5th round
    identity_rotateleft24_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [zz, LeftS_Matrices[6]]])#{6: 24}
    identity_rotateright16xor_Matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ii, zz],
        [RightS_Matrices[7], ii]])#{7: 16}
    
    
    
    def bitvectors_to_gf2vector(x, y):
        return sage.all.vector(bpr, list(int2vector(x, ws)) + list(int2vector(y, ws)))

    round_keys = [None for _ in range(51)]
    RK = get_round_keys(master_key)
    
    for i in range(len(round_keys)):
        round_keys[i] = bitvectors_to_gf2vector(RK[i].l, RK[i].r)
        #print(round_keys[i])
    
    # should be untouched Start----
    implicit_pmodadd = get_implicit_modadd_anf(ws, permuted=True, only_x_names=only_x_names)
    bpr_pmodadd = implicit_pmodadd[0].parent()
    bpr_pmodadd = BooleanPolynomialRing(names=bpr_pmodadd.variable_names(), order="deglex")
    implicit_pmodadd = [bpr_pmodadd(str(f)) for f in implicit_pmodadd]
    # should be untouched End----
    

    implicit_round_functions = []
    explicit_affine_layers = []
    i = 0
    
    affine = compose_affine(identity_matrix(2*ws),0,identity_matrix(2*ws),0)
    round_limit = 10
    for round in range( round_limit): 
        if round !=  round_limit-1:
            #print(str(i)+"\n")
            affine = compose_affine(identity_matrix(2*ws), round_keys[i],affine[0],affine[1])
            affine = compose_affine(identity_rotateright31_matrix,0,affine[0],affine[1])

            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append(affine)
            i=i+1
        
            affine = compose_affine(identity_rotateright24xor_Matrix,0,identity_rotateleft31_matrix,0)
            affine = compose_affine(identity_rotateright17_matrix,0,affine[0],0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append(affine)
            i=i+1
        
            affine = compose_affine(identity_rotateright17xor_Matrix,0,identity_rotateleft17_matrix,0)
            affine = compose_affine(identity_rotateright0_matrix,0,affine[0],0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append(affine)
            i=i+1
        
            affine = compose_affine(identity_rotateright31xor_Matrix,0,identity_rotateleft0_matrix,0)
            affine = compose_affine(identity_rotateright24_matrix,0,affine[0],0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append(affine)
            i=i+1
        
            affine = compose_affine(identity_rotateright16xor_Matrix,0,identity_rotateleft24_matrix,0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            # matrix = sage.all.block_matrix(bpr, 2, 2, [
            #         [affine[0], zero_matrix(2*ws, 2*ws)],
            #         [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            # cta = list(affine[1]) + [0 for _ in range(2*ws)]
            # anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            # if not return_implicit_round_functions:
            #     implicit_round_functions.append(anf)
            # else:
            #     implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            # if return_also_explicit_affine_layers:
            #     explicit_affine_layers.append(affine)
            i=i+1
        else:
                #print(str(i)+"\n")
            affine = compose_affine(identity_matrix(2*ws), round_keys[i],affine[0],affine[1])
            affine = compose_affine(identity_rotateright31_matrix,0,affine[0],affine[1])

            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append(affine)
            i=i+1
        
            affine = compose_affine(identity_rotateright24xor_Matrix,0,identity_rotateleft31_matrix,0)
            affine = compose_affine(identity_rotateright17_matrix,0,affine[0],0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append(affine)
            i=i+1
        
            affine = compose_affine(identity_rotateright17xor_Matrix,0,identity_rotateleft17_matrix,0)
            affine = compose_affine(identity_rotateright0_matrix,0,affine[0],0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append(affine)
            i=i+1
        
            affine = compose_affine(identity_rotateright31xor_Matrix,0,identity_rotateleft0_matrix,0)
            affine = compose_affine(identity_rotateright24_matrix,0,affine[0],0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [affine[0], zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), identity_matrix(2*ws)]])
            cta = list(affine[1]) + [0 for _ in range(2*ws)]
            anf1 = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            if return_also_explicit_affine_layers:
                explicit_affine_layers.append([affine])
            i=i+1
        
            affine = compose_affine(identity_rotateright16xor_Matrix,0,identity_rotateleft24_matrix,0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i],affine[0],0)
            affine = compose_affine(identity_matrix(2*ws),round_keys[i+1],affine[0],affine[1])
            #print(i)
            #affine = compose_affine(identity_matrix(2*ws),round_keys[i+1],affine[0],affine[1])
            aux = affine[0] ** (-1)
            matrix = sage.all.block_matrix(bpr, 2, 2, [
                    [identity_matrix(2*ws), zero_matrix(2*ws, 2*ws)],
                    [zero_matrix(2*ws, 2*ws), aux]])
            #cta = list(affine[1]) + [0 for _ in range(2*ws)]
            cta = [0 for _ in range(2*ws)] + list(aux * affine[1])
            
            anf2 = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd, bin_vector=cta)
            anf = compose_anf_fast(anf1, anf2)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
            if return_also_explicit_affine_layers:
                explicit_affine_layers[-1].append(affine)
            i=i+1
            
    print(i)        
       
        
        
        
   
    if return_also_explicit_affine_layers:
        return implicit_round_functions, explicit_affine_layers
    else:
        return implicit_round_functions


def bitvectors_to_gf2vector(x, y, ws):
    return sage.all.vector(sage.all.GF(2), list(int2vector(x, ws)) + list(int2vector(y, ws)))


def gf2vector_to_bitvectors(v, ws):
    return vector2int(v[:ws]), vector2int(v[ws:])



if __name__ == '__main__':
    
    parser = ArgumentParser(prog="sage -python crax.py", description="Generate the implicit (unencoded) and explicit affine layers of a CRAX instance for a fixed key")
    parser.add_argument("--key", nargs="+", help="the master key given as a hexadecimal representation of the words")
    parser.add_argument("--block-size", nargs="?", type=int, choices=[8, 32, 64, 128], help="the block size in bits of the CRAX instance")
    parser.add_argument("--output-file", nargs="?", help="the file to store the implicit (unencoded) and explicit affine layers")

    args = parser.parse_args()

    assert not os.path.isfile(args.output_file), f"{args.output_file} already exists"

    assert len(args.key) == 4, "key should be 4 words"
    master_key = tuple(map(lambda k: int(k, 16), args.key))
    round_keys = [None for _ in range(51)]
    bpr = sage.all.GF(2)
    rk = get_round_keys(master_key )
    def bitvectors_to_gf2vector(x, y):
        return sage.all.vector(bpr, list(int2vector(x, 32)) + list(int2vector(y, 32)))

    # i = 0
    # for K in rk:
    #     print ("round ["+str(i)+"]")
    #     print(str(hex(K.l))+" "+str(hex(K.r)))
    #     i=i+1
    
    # for i in range(len(round_keys)):
    #     round_keys[i] = bitvectors_to_gf2vector(rk[i].l, rk[i].r)
    #     print(i)
    #     print(round_keys[i])
        
    
    implicit_affine_layers, explicit_affine_layers = get_implicit_unencoded_affine_layers(master_key, return_also_explicit_affine_layers=True)
    for i, affine_layer in enumerate(implicit_affine_layers):
        # Wrap in tuple because BooleanPolynomialVector can't be pickled.
        implicit_affine_layers[i] = tuple(affine_layer)

    sage.all.save((implicit_affine_layers, explicit_affine_layers), args.output_file, compress=True)
