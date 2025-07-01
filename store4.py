import re
import numpy as np
import itertools

def monomial_index(monomial_vars, var_num, degree_offset):
    """ 计算多项式在指定阶次下的索引（用于C(n,k)展开的顺序） """
    # monomial_vars: sorted list of variable indices
    # degree_offset: 起始位置偏移，比如degree=3时从 C(n,2)+0 开始
    return degree_offset + comb_index(monomial_vars)

def comb_index(indices):
    """ 给定有序组合 [i,j,...]，返回其在C(n,k)中的字典序编号 """
    k = len(indices)
    total = 0
    last = -1
    for i, v in enumerate(indices):
        for j in range(last+1, v):
            total += comb(j, k - i - 1)
        last = v
    return total

def comb(n, k):
    from math import comb as c
    return c(n, k)

def parse_anf_to_binary_deg4(anf_str, var_num=64):
    expr = anf_str.replace(" ", "").split('+')
    
    # 计算偏移位置
    offset_deg2 = 0
    offset_deg3 = offset_deg2 + comb(var_num, 2)
    offset_deg4 = offset_deg3 + comb(var_num, 3)
    offset_deg1 = offset_deg4 + comb(var_num, 4)
    offset_const = offset_deg1 + var_num
    total_len = offset_const + 1

    bits = np.zeros(total_len, dtype=np.uint8)

    var_pattern = re.compile(r'x(\d+)')

    for term in expr:
        if term == '1':
            bits[offset_const] = 1
        else:
            matches = sorted(map(int, var_pattern.findall(term)))
            deg = len(matches)

            if deg == 1:
                bits[offset_deg1 + matches[0]] = 1
            elif deg == 2:
                bits[offset_deg2 + comb_index(matches)] = 1
            elif deg == 3:
                bits[offset_deg3 + comb_index(matches)] = 1
            elif deg == 4:
                bits[offset_deg4 + comb_index(matches)] = 1
            else:
                raise ValueError(f"不支持的阶次（{deg}次项）: {term}")

    return bits

def save_deg4_anf_vector_w(bit_vector, base_filename="anf_deg4"):
    import numpy as np
    import os

    # 保存为TXT（字符串格式）
    with open(f"{base_filename}.txt", "w") as f:
        f.write(''.join(map(str, bit_vector)))

    # 保存为 NPY（numpy 专用格式）
    np.save(f"{base_filename}.npy", bit_vector)

    # 保存为原始BIN（每 bit 用 1 字节）
    bit_vector.astype(np.uint8).tofile(f"{base_filename}_raw.bin")

    # 保存为 Bit-Packed BIN（推荐）
    pad_len = (-len(bit_vector)) % 8
    if pad_len:
        bit_vector = np.concatenate([bit_vector, np.zeros(pad_len, dtype=np.uint8)])
    packed_bytes = np.packbits(bit_vector)
    with open(f"{base_filename}_packed.bin", "wb") as f:
        f.write(packed_bytes)

    print(f"✅ 保存完成，共 {len(bit_vector)} bits → {len(bit_vector)//8} 字节（打包后）")

def save_deg4_anf_vector_a(bit_vector, base_filename="anf_deg4"):
    import numpy as np
    import os

    # 保存为TXT（字符串格式）
    with open(f"{base_filename}.txt", "a") as f:
        f.write(''.join(map(str, bit_vector)))

    # 保存为 NPY（numpy 专用格式）
    np.save(f"{base_filename}.npy", bit_vector)

    # 保存为原始BIN（每 bit 用 1 字节）
    bit_vector.astype(np.uint8).tofile(f"{base_filename}_raw.bin")

    # 保存为 Bit-Packed BIN（推荐）
    pad_len = (-len(bit_vector)) % 8
    if pad_len:
        bit_vector = np.concatenate([bit_vector, np.zeros(pad_len, dtype=np.uint8)])
    packed_bytes = np.packbits(bit_vector)
    with open(f"{base_filename}_packed.bin", "ab") as f:
        f.write(packed_bytes)

    print(f"✅ 保存完成，共 {len(bit_vector)} bits → {len(bit_vector)//8} 字节（打包后）")

# anf = "x1*x2 + x0*x1*x3 + x2*x3*x4*x5 + x7 + 1"
# bit_vector = parse_anf_to_binary_deg4(anf)
# print(f"向量长度: {len(bit_vector)}")  # 应该是 678121
# print(f"非零项数量: {np.sum(bit_vector)}")  # 应该是 5
