
import sys
import numpy as np
import binascii
import struct




file_path = 'SVDMatrix_cluster2.bin'

with open(file_path, mode='rb') as file: # b is important -> binary
    fileContent = file.read()

struct.unpack("iiiii", fileContent[:20])