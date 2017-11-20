"""

Venn test

Ideally this would be turned into a test case...

"""

from glbase3 import *

A = genelist(name="A")
B = genelist(name="B")
C = genelist(name="C")
D = genelist(name="D")
E = genelist(name="E")

A.load_list([{"key": i} for i in [1, 2, 3, 4, 5,    7, 8, 9, 10, 11, 12,                 17, 18, 19, 20,                         27]])
B.load_list([{"key": i} for i in [1, 2, 3, 4,    6, 7, 8, 9,             13, 14, 15,     17,             21, 22, 23,                 28]])
C.load_list([{"key": i} for i in [1, 2, 3,    5, 6, 7,       10, 11,     13, 14,     16,     18,         21,         24, 25,             29]])
C.load_list([{"key": i} for i in [1, 2, 3,    5, 6, 7,       10, 11,     13, 14,     16,     18,         21,         24, 25,             29]])

D.load_list([{"key": i} for i in [1, 2,    4, 5, 6,    8,    10,     12, 13,     15, 16,         19,         22,     24,     26,             30]])
E.load_list([{"key": i} for i in [1,    3, 4, 5, 6,       9,     11, 12,     14, 15, 16,             20,         23,     25, 26,                31]])       

gl = glglob(A, B, C, D) # pen and paper tested as correct, should be all 2's in each box
gl.venn(filename="test_images/test_venn4_ABCD.png", key="key")

altC = C.deepcopy()
altC.name = "altC"
gl = glglob(A, B, C, altC) 
gl.venn(filename="test_images/test_venn4_ABCC.png", key="key")

C3 = C.deepcopy()
C3.name = "C3"
gl = glglob(A, C3, C, altC) 
gl.venn(filename="test_images/test_venn4_ACCC.png", key="key")

#gl = glglob(A, B, C, D, E) # should be all 1's in each box
#gl.venn(filename="test_venn5.png", key="key")

