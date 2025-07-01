import sage.all
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

# Example: construct 128×128 block matrix using the first affine matrix
# from sage.all import BooleanPolynomialRing

bpr =sage.all.BooleanPolynomialRing(64, 'x')

affine_list = read_affine_transforms_sage_hex("build/AF.txt", bpr)

mat, vec = affine_list[0]  # use as needed
print(affine_list[0])