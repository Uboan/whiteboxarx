def ROT(x, n):
    return ((x) >> (n)) | ((x) << (32-(n)))

def ALZETTE(x, y, c):
    x += ROT(y, 31)
    y ^= ROT(x, 24)
    x ^= c
    x += ROT(y, 17)
    y ^= ROT(x, 17)
    x ^= c
    x += y
    y ^= ROT(x, 31)
    x ^= c
    x += ROT(y, 24)
    y ^= ROT(x, 16)
    x ^= c

def ALZETTE_INV(x, y, c):
    x ^= c
    y ^= ROT(x, 16)
    x -= ROT(y, 24)
    x ^= c
    y ^= ROT(x, 31)
    x -= y
    x ^= c
    y ^= ROT(x, 17)
    x -= ROT(y, 17)
    x ^= c
    y ^= ROT(x, 24)
    x -= ROT(y, 31)

NSTEPS = 1

RCON = [0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738, 0xBB1185EB]

def craxs10_enc_ref(xword, yword, key):
    for step in range(NSTEPS):
        xword ^= step
        if step % 2 == 0:
            xword ^= key[0]
            yword ^= key[1]
        else:
            xword ^= key[2]
            yword ^= key[3]
        ALZETTE(xword[0], yword[0], RCON[step % 5])
    xword[0] ^= key[0]
    yword[0] ^= key[1]

def craxs10_dec_ref(xword, yword, key):
    xword[0] ^= key[0]
    yword[0] ^= key[1]
    for step in range(NSTEPS-1, -1, -1):
        ALZETTE_INV(xword[0], yword[0], RCON[step % 5])
        if step % 2 == 0:
            xword[0] ^= key[0]
            yword[0] ^= key[1]
        else:
            xword[0] ^= key[2]
            yword[0] ^= key[3]
        xword[0] ^= step

# Example usage
plaintext = [0x0f0e0d0c, 0x0b0a0908]
key = [0x0f0e0d0c, 0x0b0a0908, 0x07060504, 0x03020100]

# Encryption
craxs10_enc_ref(plaintext, plaintext[1:], key)

# Decryption
#craxs10_dec_ref(plaintext, plaintext[1:], key)

# Print the result
print("Plaintext after encryption and decryption:", plaintext)
