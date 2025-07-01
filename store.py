import re
import numpy as np

def parse_anf_to_binary(anf_str, var_num=64, repeat=2):
    expr = anf_str.replace(" ", "").split('+')

    # 二次项的映射：x0x1, x0x2, ..., x62x63（共2016个）
    quad_map = {}
    idx = 0
    for i in range(var_num):
        for j in range(i+1, var_num):
            quad_map[(i, j)] = idx
            idx += 1
    num_quad = idx  # =2016
    num_linear = var_num
    total_len = num_quad + num_linear + 1

    # 初始化结果向量
    bits = np.zeros(total_len, dtype=np.uint8)

    # 正则用于匹配 x<num>
    var_pattern = re.compile(r'x(\d+)')

    for term in expr:
        if term == '1':
            bits[-1] = 1  # 常数项
        else:
            # 提取所有变量编号
            matches = var_pattern.findall(term)
            if len(matches) == 2:
                i, j = sorted(map(int, matches))
                bits[quad_map[(i, j)]] = 1
            elif len(matches) == 1:
                i = int(matches[0])
                bits[num_quad + i] = 1
            else:
                raise ValueError(f"无法识别的项: {term}")

    # 重复2次输出
    return np.tile(bits, repeat)

def save_as_bitpacked_bin_w(bit_vector, filename):#覆盖文件内容
    # 确保是 0/1 的 uint8 向量
    bit_vector = np.array(bit_vector, dtype=np.uint8)

    # 如果长度不是8的倍数，填充0
    pad_len = (-len(bit_vector)) % 8
    if pad_len:
        bit_vector = np.concatenate([bit_vector, np.zeros(pad_len, dtype=np.uint8)])

    # 将每8位合并成一个字节
    packed_bytes = np.packbits(bit_vector)

    # 写入文件
    with open(filename, "wb") as f:#覆盖文件内容
        f.write(packed_bytes)
        
def save_as_bitpacked_bin_a(bit_vector, filename):#在文件后面追加
    import os
    # 确保是 0/1 的 uint8 向量
    bit_vector = np.array(bit_vector, dtype=np.uint8)

    # 如果长度不是8的倍数，填充0
    pad_len = (-len(bit_vector)) % 8
    if pad_len:
        bit_vector = np.concatenate([bit_vector, np.zeros(pad_len, dtype=np.uint8)])

    # 将每8位合并成一个字节
    packed_bytes = np.packbits(bit_vector)
    # print("+++++++++++packed_bytes+++++++++++++")
    # print(len(packed_bytes))
    # print("+++++++++++packed_bytes+++++++++++++")
    # 写入文件
    with open(filename, "ab") as f:#在文件后面追加
        f.write(packed_bytes)
        
    # print("+++++++++++packed_bytes+++++++++++++")
    # print(len(packed_bytes))
    # print("+++++++++++packed_bytes+++++++++++++")
    print("+++++++++++文件：{filename} 字节+++++++++++++")
    size = os.path.getsize(filename)
    print(f"文件大小：{size} 字节")
    print("+++++++++++文件大小：size 字节+++++++++++++")
    return packed_bytes

# 示例字符串
# anf_expr = """x0*x3 + x1*x2 + x4*x5 + x5*x6 + x6*x7 + x8*x11 + x12*x13 + x12*x15 + x13*x14 + x14*x15 + 
# x16*x18 + x17*x18 + x20*x22 + x20*x23 + x21*x22 + x24*x25 + x24*x27 + x25*x26 + x26*x27 + x28*x29 + 
# x29*x30 + x30*x31 + x32*x33 + x32*x34 + x32*x35 + x34*x35 + x36*x37 + x36*x38 + x36*x39 + x37*x38 + 
# x38*x39 + x40*x42 + x41*x42 + x44*x47 + x48*x49 + x48*x50 + x48*x51 + x50*x51 + x52*x53 + x52*x54 + 
# x52*x55 + x53*x54 + x54*x55 + x56*x57 + x57*x58 + x58*x59 + x60*x62 + x0 + x1 + x3 + x4 + x8 + x9 + 
# x10 + x13 + x14 + x16 + x17 + x19 + x22 + x23 + x25 + x26 + x28 + x32 + x34 + x35 + x38 + x41 + x42 + 
# x43 + x44 + x45 + x46 + x48 + x50 + x51 + x54 + x58 + x60 + x61 + x62 + 1"""

# # 调用函数
# bit_vector = parse_anf_to_binary(anf_expr)

# # 保存为文本（推荐阅读）
# with open("anf_bits.txt", "w") as f:
#     f.write(''.join(map(str, bit_vector)))

# # 保存为npy（推荐用于Python分析）
# np.save("anf_bits.npy", bit_vector)

# # 保存为真正的二进制文件（推荐用于嵌入式/底层）
# bit_vector.astype(np.uint8).tofile("anf_bits.bin")
# save_as_bitpacked_bin(bit_vector, "anf_bits_packed.bin")
