"""
expression.py tester code.

Part of glbase

Tests for expression class

"""

import unittest

# get glbase
import sys, os, math
sys.path.append(os.path.realpath("../../"))

import glbase as gl

def fq_eq(a, b, eps=0.000001):
    return abs(math.log(abs(a)) - math.log(abs(b))) <=  eps

class Test_Expression(unittest.TestCase):
    def setUp(self):
        data = [{"name": "i1", "conditions": [1, 2,    3, 4], "m1": "aaa", "m2": "bb"},
                {"name": "i2", "conditions": [1, 1,    2, 1], "m1": "zzz", "m2": "aa"},
                {"name": "i3", "conditions": [4, 3,    2, 1], "m1": "aaa", "m2": "aa"},
                {"name": "i4", "conditions": [1, 2.05, 3, 0.1], "m1": "zzz", "m2": "cc"}]
        
        data_err = [{"name": "i1", "conditions": [1, 2,    3, 4], "err": [1, 2,    3, 4], "m1": "aaa", "m2": "bb"},
            {"name": "i2", "conditions": [1, 1,    1, 1], "err": [1, 1,    1, 1], "m1": "zzz", "m2": "aa"},
            {"name": "i3", "conditions": [4, 3,    2, 1], "err": [4, 3,    2, 1], "m1": "aaa", "m2": "aa"},
            {"name": "i4", "conditions": [1, 2.05, 3, 0.1], "err": [1, 2.05, 3, 0.1], "m1": "zzz", "m2": "cc"}]
        
        self.expn = gl.expression(loadable_list=data, cond_names=["a", "b", "c", "d"])
        self.expn_err = gl.expression(loadable_list=data_err, cond_names=["a", "b", "c", "d"])

    def test_fold_change(self):
        
        expn = self.expn.add_fc_key(key="fca", cond1="a", cond2="b")
        self.assertEqual(expn[0]["fca"], 0.9999992786530209) # not 2.0 because of pad
        self.assertEqual(expn[1]["fca"], -0.0)
        self.assertEqual(expn[2]["fca"], -0.41503737905429205)
        self.assertEqual(expn[3]["fca"], 1.0356231707899088)
        
        tt = expn.norm_multi_fc({"a": ["b"], "c": ["d"]}, log=None) 
        # These test a vs b
        self.assertEqual(tt[0]["conditions"][1], 1.9999990000010004) # not 2.0 because of pad
        self.assertEqual(tt[1]["conditions"][1], -1.0)
        self.assertEqual(tt[2]["conditions"][1], -1.3333332222222591)
        self.assertEqual(tt[3]["conditions"][1], 2.0499989500010503)
        
        self.assertEqual(tt[0]["conditions"][0], 0.0) # test blank gets loaded.
        
        # test c vs d is distinct from a vs d:
        self.assertEqual(tt[0]["conditions"][3], 1.3333332222222591) #
        self.assertEqual(tt[1]["conditions"][3], -1.9999990000010004)
        self.assertEqual(tt[2]["conditions"][3], -1.9999990000010004)
        self.assertEqual(tt[3]["conditions"][3], -29.99971000289997)   
                
        self.assertEqual(tt[0]["conditions"][0], 0.0) # test blank gets loaded.     

    def test_multi_sort(self):
        # multi_sort is a bit tricky in expression objects. This is a test case.
        self.expn.multi_sort(["m1", "m2"]) # keys only
        exp_order = ["i3", "i1", "i2", "i4"]
        obs_order = [i["name"] for i in self.expn] 
        self.assertTrue(False not in [i1 == i2 for i1, i2 in zip(obs_order, exp_order)])

        self.expn.multi_sort(["a", "b"]) # conditions only      
        exp_order = ["i2", "i1", "i4", "i3"]
        obs_order = [i["name"] for i in self.expn] 
        self.assertTrue(False not in [i1 == i2 for i1, i2 in zip(obs_order, exp_order)])
        
        self.expn.multi_sort(["m1", "c"]) # keys and conditions
        exp_order = ["i3", "i1", "i2", "i4"]
        obs_order = [i["name"] for i in self.expn] 
        self.assertTrue(False not in [i1 == i2 for i1, i2 in zip(obs_order, exp_order)])

    def test_slice_conditions(self):
        slic = self.expn.sliceConditions(["a", "c"])
        self.assertEqual(len(slic.getConditionNames()), 2)
        self.assertEqual(len(slic), 4)
        self.assertEqual(slic[0]["conditions"], [1.0, 3.0])
        self.assertEqual(slic[2]["conditions"], [4.0, 2.0])

        slic = self.expn.sliceConditions(["b", "d", "c"])
        self.assertEqual(len(slic.getConditionNames()), 3)
        self.assertEqual(len(slic), 4)
        self.assertEqual(slic[0]["conditions"], [2.0, 4.0, 3.0])
        self.assertEqual(slic[2]["conditions"], [3.0, 1.0, 2.0])

    def test_slice_conditions_errors(self):
        slic = self.expn_err.sliceConditions(["a", "c"])
        self.assertEqual(len(slic.getConditionNames()), 2)
        self.assertEqual(len(slic), 4)
        self.assertEqual(slic[0]["err"], [1.0, 3.0])
        self.assertEqual(slic[2]["err"], [4.0, 2.0])

        slic = self.expn_err.sliceConditions(["b", "d", "c"])
        self.assertEqual(len(slic.getConditionNames()), 3)
        self.assertEqual(len(slic), 4)
        self.assertEqual(slic[0]["err"], [2.0, 4.0, 3.0])
        self.assertEqual(slic[2]["err"], [3.0, 1.0, 2.0])

    def test_reloading_from_numpy_array(self):
        self.expn.subtract_mean() # Checked by hand
        self.assertEqual(self.expn[0]["conditions"], [-0.75, -0.012500000000000178, 0.5, 2.4750000000000001])
        
        self.expn.whiten() # I'm just assuming this one is correct, it's not clear how to check it.
        self.assertEqual(self.expn[0]["conditions"], [-0.57735026918962573, -0.017669388943904119, 1.0, 1.6774842736586513])
    
    def test_normalization_cols(self):
        self.expn.normalize_columns()
        self.assertEqual(self.expn[1]["conditions"], [-1.0, -1.0, -1.0, -0.53846153846153855])       
        
    #def test_normalization_rows(self): 
    #    print self.expn.all()
    #    self.expn.normalize_rows()
    #    print self.expn.all()

    def test_digitize(self): 
        self.expn.digitize(8)
        self.assertEqual(self.expn[0]["conditions"], [2.0, 4.0, 6.0, 8.0])       
        #self.expn.digitize(4)
        self.assertEqual(self.expn[1]["conditions"], [2.0, 2.0, 4.0, 2.0])   
    
    def test_stats(self):
        # test fq_eq is appropriate for this tollerances
        self.assertFalse(fq_eq(0.001, 0.00011))
        self.assertTrue(fq_eq(1.000001, 1.0000012))
        t = self.expn.stats.ttest('a', 'b')
        self.assertTrue(fq_eq(t[0], -0.3073755916437676))
        self.assertTrue(fq_eq(t[1], 0.77188223882169515))
        m = self.expn.stats.min()
        self.assertTrue(fq_eq(m, 0.1))
        m = self.expn.stats.max()
        self.assertTrue(fq_eq(m, 4.0))
        a = self.expn.stats.mean()
        self.assertTrue(fq_eq(a, 1.9468750000000001)) 
        s = self.expn.stats.stdev()
        self.assertTrue(fq_eq(s, 1.1337808361297168))
        p = self.expn.stats.pearsonr('a', 'b')
        self.assertTrue(fq_eq(p[0], 0.80591269059114867))
        self.assertTrue(fq_eq(p[1], 0.19408730940885136))
        #p = self.expn.stats.print_stats()
        
    def test_map_with_homogenous_and_heterogenous_condition_names(self):
        A = [{"name": "gene1", "conditions": [1,2,3,4]}]
        B = [{"name": "gene1", "conditions": [7,8,9,0]}]

        expn_A = gl.expression(loadable_list=A, cond_names=["cA", "cB", "cC", "cD"])
        expn_B = gl.expression(loadable_list=B, cond_names=["cA", "cB", "cC", "cD"]) # same copndition names

        m = expn_A.map(genelist=expn_B, key="name")
        
        self.assertEqual(m[0]["name"], "gene1")
        self.assertEqual(m[0]["conditions"], [7.0, 8.0, 9.0, 0.0])
        self.assertEqual(m.getConditionNames(), ['cA', 'cB', 'cC', 'cD'])

        A = [{"name": "gene1", "conditions": [1,2,3,4]}]
        C = [{"name": "gene1", "conditions": [7,8,9,0]}]

        expn_A = gl.expression(loadable_list=A, cond_names=["cA", "cB", "cC", "cD"])
        expn_C = gl.expression(loadable_list=B, cond_names=["cE", "cF", "cG", "cH"]) # now we have differnt names

        m = expn_A.map(genelist=expn_C, key="name")

        self.assertEqual(m[0]["name"], "gene1")
        self.assertEqual(m[0]["conditions"], [7.0, 8.0, 9.0, 0.0, 1.0, 2.0, 3.0, 4.0])
        self.assertEqual(m.getConditionNames(), ['cE', 'cF', 'cG', 'cH', 'cA', 'cB', 'cC', 'cD'])

    def test_correlation_heatmap(self):
        data = self.expn.correlation_heatmap(axis="conditions", filename="test_images/corr_heatmap1.png", mode="r2")
        data = self.expn.correlation_heatmap(axis="conditions", filename="test_images/corr_heatmap2.png", mode="pearson")
        data = self.expn.correlation_heatmap(axis="conditions", filename="test_images/corr_heatmap3.png", mode="spearman")
        data = self.expn.correlation_heatmap(axis="genes", label_key="name", filename="test_images/corr_heatmap4.png", mode="r2")
        data = self.expn.correlation_heatmap(axis="genes", label_key="name", filename="test_images/corr_heatmap5.png", mode="pearson")
        data = self.expn.correlation_heatmap(axis="genes", label_key="name", filename="test_images/corr_heatmap6.png", mode="spearman")
    
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Expression)
    unittest.TextTestRunner(verbosity=2).run(suite)
