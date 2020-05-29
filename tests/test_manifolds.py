"""
manifolds tester code.

Part of glbase3

Tests for the various manifolds

"""

import unittest

# get glbase
import sys, os, math
sys.path.append(os.path.realpath("../../"))

import glbase3 as gl

def fq_eq(a, b, eps=0.000001):
    return abs(math.log(abs(a)) - math.log(abs(b))) <=  eps

class Test_Manifold(unittest.TestCase):
    def setUp(self):
        data = [{"name": "i1",  "conditions": [1, 2, 3, 4]},
                {"name": "i2",  "conditions": [1, 1, 2, 1]},
                {"name": "i3",  "conditions": [4, 3, 2, 1]},
                {"name": "i4",  "conditions": [1, 2, 3, 1]},
                {"name": "i5",  "conditions": [2, 1, 5, 1]},
                {"name": "i6",  "conditions": [3, 3, 3, 1]},
                {"name": "i7",  "conditions": [4, 2, 3, 1]},
                {"name": "i8",  "conditions": [2, 1, 5, 1]},
                {"name": "i9",  "conditions": [3, 2, 3, 1]},
                {"name": "i10", "conditions": [4, 3, 3, 1]},
                {"name": "i11", "conditions": [1, 1, 5, 1]},
                ]

        self.expn = gl.expression(loadable_list=data, cond_names=["a", "b", "c", "d"])

    def test_tsne(self):
        tsne = self.expn.tsne.configure(whiten=True, random_state=42, verbose=False)
        self.assertTrue(self.expn.tsne.whiten)
        self.expn.tsne.train(2)
        self.assertTrue(self.expn.tsne.trained)
        self.expn.tsne.scatter(filename='/tmp/tsne_scat.png')
        ret = self.expn.tsne.cluster(num_clusters=2, method='KMeans', filename='/tmp/tsne_scat_clus.png')

        self.assertListEqual(list(ret[1]), [1, 0, 0, 1])
        self.assertTrue(fq_eq(ret[2][0][0], -73.73201))
        self.assertTrue(fq_eq(ret[2][1][1], -51.5523 ))

    def test_mds(self):
        mds = self.expn.mds.configure(whiten=True, random_state=42, verbose=False)
        self.assertTrue(self.expn.mds.whiten)
        self.expn.mds.train(2)
        self.assertTrue(self.expn.mds.trained)
        self.expn.mds.scatter(filename='/tmp/mds_scat.png')
        ret = self.expn.mds.cluster(num_clusters=2, method='KMeans', filename='/tmp/mds_scat_clus.png')

        print(ret)
        self.assertListEqual(list(ret[1]), [1, 1, 0, 1])
        self.assertTrue(fq_eq(ret[2][0][0], -0.52325895))
        self.assertTrue(fq_eq(ret[2][1][1], 0.449152 ))

    def test_umap(self):
        if gl.config.UMAP_LEARN_AVAIL:
            mds = self.expn.mds.configure(whiten=True, random_state=42, verbose=False)
            self.assertTrue(self.expn.mds.whiten)
            self.expn.mds.train(2)
            self.assertTrue(self.expn.mds.trained)
            self.expn.mds.scatter(filename='/tmp/mds_scat.png')
            ret = self.expn.mds.cluster(num_clusters=2, method='KMeans', filename='/tmp/mds_scat_clus.png')

            print(ret)
            self.assertListEqual(list(ret[1]), [1, 1, 0, 1])
            self.assertTrue(fq_eq(ret[2][0][0], -0.52325895))
            self.assertTrue(fq_eq(ret[2][1][1], 0.449152 ))

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Manifold)
    unittest.TextTestRunner(verbosity=2).run(suite)
