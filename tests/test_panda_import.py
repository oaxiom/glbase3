

import glbase3 as gl
import pandas as pd
import numpy as np

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase3

class Test_PandasImportlist(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame(
            {'A' : 1.,
            'C' : pd.Series(1,index=list(range(4)),dtype='float32'),
            'D' : [1, 2, 3, 4],
            'E' : pd.Categorical(["test","train","test","train"]),
            'F' : 'foo'})
        
        self.dfexpn = pd.DataFrame(
            {'name': ["Sox2","Nanog","Esrrb","Stat3"],
            'ensg': ["ENSG1","ENSG2","ENSG3","ENSG4"],
            'expn1': [1, 3, 4, 5],
            'expn2': [4, 3, 2, 1],
            'expn3': [3, 3, 3, 3],
            })
    
    def test_genelist_from_pandas(self):
        glb = gl.genelist()
        glb.from_pandas(self.df)
       
        self.assertEqual(len(glb), self.df.shape[0])
        self.assertSetEqual(set(glb[0].values()), set(self.df.loc[[0]].values[0]))
        self.assertSetEqual(set(glb[1].values()), set(self.df.loc[[1]].values[0]))
        self.assertSetEqual(set(glb[3].values()), set(self.df.loc[[3]].values[0]))

    def test_expression_from_pandas(self):
        expn = gl.expression()
        expn.from_pandas(self.dfexpn, condition_names=['expn1', 'expn2', 'expn3'])
       
        self.assertEqual(len(expn), self.dfexpn.shape[0])
        
        self.assertListEqual(list(self.dfexpn['expn1'].values), list(expn['expn1']))
        self.assertListEqual(list(self.dfexpn['expn2'].values), list(expn['expn2']))
        self.assertListEqual(list(self.dfexpn['name'].values), list(expn['name']))
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_PandasImportlist)
    unittest.TextTestRunner(verbosity=2).run(suite)