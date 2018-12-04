import numpy as np
from sys import argv
f1 = file(argv[1])
f2 = file(argv[2])
#print f1,f2
while True:
    line1 = f1.readline()
    line2 = f2.readline()
    if len(line1) == 0 or len(line2) == 0:
        break
    words1 = line1.split()
    words2 = line2.split()
    print words1[0], \
          float(words1[1]) - float(words2[1]), \
          float(words1[2]) - float(words2[2]), \
          float(words1[3]) - float(words2[3])
