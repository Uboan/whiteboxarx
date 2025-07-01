
from collections import namedtuple
from functools import partial

import os
import sys
import sage.all
from store import * # degree 2的存储代码
from store4 import * # degree 4的存储代码
from pympler import asizeof


from boolcrypt.utilities import (
    substitute_variables, BooleanPolynomialRing, vector2int,
    int2vector, compose_affine, matrix2anf, compose_anf_fast, get_smart_print,get_ct_coeff,get_symbolic_coeff
)

from boolcrypt.modularaddition import get_4bits_different_sbox_anf

from argparse import ArgumentParser

class Present:
    default_rounds = 31
    # Shift=[31,24,17,17,0,31,24,16]
    state_size = 64
    
class Pair:
    def __init__(self):
        self.l = 0
        self.r = 0
    def __init__(self, left, right):
        self.l = left
        self.r = right
 
# P Box permutations table   
pBox = [
    0, 16, 32, 48, 1, 17, 33, 49,
    2, 18, 34, 50, 3, 19, 35, 51,
    4, 20, 36, 52, 5, 21, 37, 53,
    6, 22, 38, 54, 7, 23, 39, 55,
    8, 24, 40, 56, 9, 25, 41, 57,
    10, 26, 42, 58, 11, 27, 43, 59,
    12, 28, 44, 60, 13, 29, 45, 61,
    14, 30, 46, 62, 15, 31, 47, 63
]



def key_initial():
    RK = []
    RK.append(0x0000000000000000)
    RK.append(0xC000000000000000)
    RK.append(0x5000180000000001)
    RK.append(0x60000A0003000001)
    RK.append(0xB0000C0001400062)
    RK.append(0x900016000180002A)
    RK.append(0x0001920002C00033)
    RK.append(0xA000A0003240005B)
    RK.append(0xD000D4001400064C)
    RK.append(0x30017A001A800284)
    RK.append(0xE01926002F400355)
    RK.append(0xF00A1C0324C005ED)
    RK.append(0x800D5E014380649E)
    RK.append(0x4017B001ABC02876)
    RK.append(0x71926802F600357F)
    RK.append(0x10A1CE324D005EC7)
    RK.append(0x20D5E21439C649A8)
    RK.append(0xC17B041ABC428730)
    RK.append(0xC926B82F60835781)
    RK.append(0x6A1CD924D705EC19)
    RK.append(0xBD5E0D439B249AEA)
    RK.append(0x07B077ABC1A8736E)
    RK.append(0x426BA0F60EF5783E)
    RK.append(0x41CDA84D741EC1D5)
    RK.append(0xF5E0E839B509AE8F)
    RK.append(0x2B075EBC1D0736AD)
    RK.append(0x86BA2560EBD783AD)
    RK.append(0x8CDAB0D744AC1D77)
    RK.append(0x1E0EB19B561AE89B)
    RK.append(0xD075C3C1D6336ACD)
    RK.append(0x8BA27A0EB8783AC9)
    RK.append(0x6DAB31744F41D700)
    return RK





def read_affine_transforms_sage_hex(filename, bpr):
    """
    Reads affine transforms from a file with 16-digit hex rows and returns a list of (matrix, vector)
    using Sage's BooleanPolynomialRing.
    
    :param filename: Path to the text file
    :param bpr: A BooleanPolynomialRing(64, ...)
    :return: List of (matrix, vector) tuples
    """
    from sage.all import Matrix, vector, Integer

    def hex_to_bitlist(hex_str):
        # Convert hex string to 64-bit binary string, then to list of ints
        return [int(b) for b in bin(Integer(hex_str, 16))[2:].zfill(64)]

    transforms = []

    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    assert len(lines) % 65 == 0, "Each affine transform should be 65 lines (64 matrix + 1 vector)"

    for i in range(0, len(lines), 65):
        block = lines[i:i+65]
        matrix_rows = [hex_to_bitlist(row) for row in block[:64]]
        affine_vector = hex_to_bitlist(block[64])

        mat = Matrix(bpr, matrix_rows)
        vec = vector(bpr, affine_vector)

        transforms.append((mat, vec))

    return transforms

def read_sboxes(filename):
    sboxes = []
    with open(filename, 'r') as f:
        for line in f:
            nums = list(map(int, line.strip().split()))
            sboxes.append(nums)
    return sboxes

def __sizeof__(self) -> int:
        total = sys.getsizeof(self.myArray)
        for i in range(self.present):
            total += sys.getsizeof(self.myArray[i])
        
        return total

def get_implicit_encoded_present_functions(
        master_key,output_file_name, only_x_names=False
):
    
    rounds = 31

    
    # generating shift operation in Matrices{Step 1:In 1 word}
    bpr = sage.all.GF(2)
    def int2bitvector(input):
        vector = [0 for i in range(Present.state_size)]
        vector[input] = 1
        return vector
    

    def bitvectors_to_gf2vector(x):
        return sage.all.vector(bpr, list(int2vector(x, Present.state_size)))

    identity_matrix = partial(sage.all.identity_matrix, bpr)
    # zero_matrix = partial(sage.all.zero_matrix, bpr)
    id_round_matrix = identity_matrix(Present.state_size)

   
    pM = []
    for i in range(64):
       v = int2bitvector(pBox[i])
       pM.append(v) 
    
    P_round_matrix = sage.all.matrix(pM)
    
    round_keys = [None for _ in range(Present.default_rounds)]
    RK = key_initial()
    
    
    for i in range(len(round_keys)):
        round_keys[i] = bitvectors_to_gf2vector(RK[i])
   
    affine_encodings = read_affine_transforms_sage_hex("build/AF.txt",bpr)
    inv_affine_encodings = read_affine_transforms_sage_hex("build/AF_inv.txt",bpr)
    
    sboxes1 = read_sboxes("build/S4_D2Sbox1.txt")
    sboxes2 = read_sboxes("build/S4_D2Sbox2.txt")
    
    
    implicit_round_functions = []
    # explicit_affine_layers = []
    
    ### round 1 start --------------
    i = 0
    implicit_psboxanf = get_4bits_different_sbox_anf(Present.state_size,sboxes1[i:i+16], only_x_names=True)
    bpr_psboxanf= implicit_psboxanf[0].parent()
    bpr_psboxanf = BooleanPolynomialRing(names=bpr_psboxanf.variable_names(), order="deglex")
    implicit_psboxanf = [bpr_psboxanf(str(f)) for f in implicit_psboxanf]
    
    # print(implicit_psboxanf)

    cta = list(round_keys[i]) + [0 for _ in range(Present.state_size)]
    anf = matrix2anf(id_round_matrix, bool_poly_ring=bpr_psboxanf, bin_vector=cta)
    
    anf = compose_anf_fast(implicit_psboxanf,anf)
    
    affine_matrix,affine_vector = affine_encodings[i]
    cta = list(affine_vector) + [0 for _ in range(Present.state_size)]
    encanf = matrix2anf(affine_matrix, bool_poly_ring=bpr_psboxanf, bin_vector=cta)

    anf = list(sage.all.vector(bpr_psboxanf,compose_anf_fast(encanf,anf)))
    implicit_round_functions.append(anf)
    ### round  1  end --------------
    
   
  
    ##存储第一轮
    stranf = str(anf[0])

    bit_v = parse_anf_to_binary_deg4(stranf)

    bit_v.astype(np.uint8).tofile(output_file_name+".bin")
    # print("------------size of bit_v------------")
    # print((bit_v.nbytes))
    # print("------------size of bit_v------------")
    packed_bits = save_as_bitpacked_bin_w(bit_v,output_file_name+"packed.bin")
    
    for i in range(1,64):
        
        stranf = str(anf[i])
        # print("-----------------------")
        # print(stranf)
        bit_v = parse_anf_to_binary_deg4(stranf)
        # print("------------bit_v------------")
        # print(len(bit_v))
        # print("------------bit_v------------")
        bit_v.astype(np.uint8).tofile(output_file_name+".bin")
        # print("------------size of bit_v------------")
        # print((bit_v.nbytes))
        # print("------------size of bit_v------------")
        
        # if i==0:
            
        #     packed_bits = save_as_bitpacked_bin_w(bit_v,output_file_name+"packed.bin")
        # else:
        packed_bits = save_as_bitpacked_bin_a(bit_v,output_file_name+"packed.bin")
        print("------------size of packed_bits------------")
        print(packed_bits.nbytes)
        print("------------size of packed_bits------------")
        
        # save_deg4_anf_vector_a(bit_v,output_file_name)
        # print("=======================")
        # stranfs.append(stranf)
    
    ##存储第一轮完毕
 
    
    
    ###开始中间轮
    for i in range(1,rounds-1):
        print("round"+str(i))
        implicit_psboxanf = get_4bits_different_sbox_anf(Present.state_size,sboxes2[(i-1)*16:(i-1)*16+16], only_x_names=only_x_names)
        bpr_psboxanf= implicit_psboxanf[0].parent()
        bpr_psboxanf = BooleanPolynomialRing(names=bpr_psboxanf.variable_names(), order="deglex")
        implicit_psboxanf = [bpr_psboxanf(str(f)) for f in implicit_psboxanf]
        print("round"+str(i)+" Sbox end")
        # compose affineU64 and S2  
        affine_matrix,affine_vector = inv_affine_encodings[i-1]
        cta = list(affine_vector) + [0 for _ in range(Present.state_size)]
        inv_encanf = matrix2anf(affine_matrix, bool_poly_ring=bpr_psboxanf, bin_vector=cta)
        anf = compose_anf_fast(implicit_psboxanf,inv_encanf)
        # anf = inv_encanf # test
        print("round"+str(i)+" affineU64 end")
        
        
        # compose it with P and AddRoundKey
        cta = list(round_keys[i]) + [0 for _ in range(Present.state_size)]
        p_anf = matrix2anf(P_round_matrix, bool_poly_ring=bpr_psboxanf, bin_vector=cta)
        anf = compose_anf_fast(p_anf,anf)
        # anf = p_anf # test
        
        print("round"+str(i)+" compose it with P and AddRoundKey end")
        # S1
        implicit_psboxanf = get_4bits_different_sbox_anf(Present.state_size,sboxes1[(i)*16:(i)*16+16], only_x_names=only_x_names)
        bpr_psboxanf= implicit_psboxanf[0].parent()
        bpr_psboxanf = BooleanPolynomialRing(names=bpr_psboxanf.variable_names(), order="deglex")
        implicit_psboxanf = [bpr_psboxanf(str(f)) for f in implicit_psboxanf]
        
        print("round"+str(i)+" S1 end")
        
        # compose the previous and S1
        anf = compose_anf_fast(implicit_psboxanf,anf)
        
        print("round"+str(i)+" compose the previous and S1 end")
        #affineU64  
        affine_matrix,affine_vector = affine_encodings[i]
        cta = list(affine_vector) + [0 for _ in range(Present.state_size)]
        encanf = matrix2anf(affine_matrix, bool_poly_ring=bpr_psboxanf, bin_vector=cta)
        print("round"+str(i)+" affineU64 end")
        # print(encanf)
        # compose affineU64 and the previous
        anf = list(sage.all.vector(bpr_psboxanf,compose_anf_fast(encanf,anf)))
        print("round"+str(i)+" compose affineU64 and the previous end")
        anf = list(sage.all.vector(bpr_psboxanf,anf))
        # print(anf)
        implicit_round_functions.append(anf)
        
        ###开始存储中间轮
        for ip in range(64):
            stranf = str(anf[ip])
            print("--round"+str(i)+"-------"+str(ip)+"--------------")
            # print(stranf)
            bit_v = parse_anf_to_binary_deg4(stranf)
            
            bit_v.astype(np.uint8).tofile(output_file_name+".bin")
            save_as_bitpacked_bin_a(bit_v,output_file_name+"packed.bin")

    
    implicit_psboxanf = get_4bits_different_sbox_anf(Present.state_size,sboxes2[(rounds-1)*16:(rounds-1)*16+16], only_x_names=only_x_names)
    bpr_psboxanf= implicit_psboxanf[0].parent()
    bpr_psboxanf = BooleanPolynomialRing(names=bpr_psboxanf.variable_names(), order="deglex")
    implicit_psboxanf = [bpr_psboxanf(str(f)) for f in implicit_psboxanf]
    
    # compose affineU64 and S2  
    affine_matrix,affine_vector = inv_affine_encodings[rounds-1]
    cta = list(affine_vector) + [0 for _ in range(Present.state_size)]
    inv_encanf = matrix2anf(affine_matrix, bool_poly_ring=bpr_psboxanf, bin_vector=cta)
    anf = compose_anf_fast(implicit_psboxanf,inv_encanf)
    
    # compose it with AddRoundKey
    cta = list(round_keys[rounds-1]) + [0 for _ in range(Present.state_size)]
    p_anf = matrix2anf(id_round_matrix, bool_poly_ring=bpr_psboxanf, bin_vector=cta)
    # implicit_round_functions.append(compose_anf_fast(p_anf,anf))
    
    anf = list(sage.all.vector(bpr_psboxanf,compose_anf_fast(encanf,anf)))
    implicit_round_functions.append(anf) 

    
    
    for i in range(64):
        stranf = str(anf[i])
        # print("-----------------------")
        # print(stranf)
        bit_v = parse_anf_to_binary_deg4(stranf)
        # save_as_bitpacked_bin_a(bit_v,output_file_name)
        
        bit_v.astype(np.uint8).tofile(output_file_name+".bin")
        save_as_bitpacked_bin_a(bit_v,output_file_name+"packed.bin")

        
    return implicit_round_functions


def bitvectors_to_gf2vector(x, y, ws):
    return sage.all.vector(sage.all.GF(2), list(int2vector(x, ws)) + list(int2vector(y, ws)))


def gf2vector_to_bitvectors(v, ws):
    return vector2int(v[:ws]), vector2int(v[ws:])



if __name__ == '__main__':
    
    parser = ArgumentParser(prog="sage -python present.py", description="Generate the implicit encoded a present instance for a fixed key")
    parser.add_argument("--output-file", nargs="?", help="the file to store the implicit encoded instance")

    args = parser.parse_args()

    assert not os.path.isfile(args.output_file), f"{args.output_file} already exists"
    bpr = sage.all.GF(2)

        
    output_file_name = "31round_anf_bits"
    implicit_functions = get_implicit_encoded_present_functions(None,output_file_name,True)
    
    
    # print("sizeof item in implicit_functions")
    # for it in implicit_functions:
    #     # Wrap in tuple because BooleanPolynomialVector can't be pickled.
    #     it = tuple(it)
    #     print(asizeof.asizeof(it))
    # print("what???????")
    # print(asizeof.asizeof(implicit_functions))
    # sage.all.save((implicit_functions), args.output_file)
