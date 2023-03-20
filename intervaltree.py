"""

This was adapted from an old py2.7, pure python bx-python:

Copyright (c) 2005-2015 The Pennsylvania State University
Copyright (c) 2013-2015 The Johns Hopkins University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

https://bitbucket.org/james_taylor/bx-python/raw/ebf9a4b352d3/lib/bx/intervals/operations/quicksect.py

"""
import math
import time
import sys
import random

class intervaltree(object):
    def __init__(self):
        """
        **Purpose**
            Entry class for the intervaltree, a fast method of comparing and intersecting
            intervals

        **Returns**
            None
        """
        self.chroms = {}
        self.chroms_list = []

    def insert(self, chrom, left, right, linenum=0, other=None):
        if chrom in self.chroms:
            self.chroms[chrom] = self.chroms[chrom].insert(left, right, linenum, other)
        else:
            self.chroms[chrom] = intervalnode(left, right, linenum, other)
        self.chroms_list = [(start, end, False) for start, end in self.chroms_list[chrom]]

    def intersect( self, interval, report_func ):
        chrom = interval.chrom
        if chrom in self.chroms:
            start = interval.start
            end = interval.end
            self.chroms[chrom].intersect( start, end, report_func )

    def traverse( self, func ):
        for item in self.chroms.itervalues():
            item.traverse( func )

    def maximize_nonoverlapping_count(self, chrom):
        '''
        **Purpose**
            For a list of intervals, work out the maximal nonoverlapping set

        **Arguments**
            chrom (Required)
                Chromosome name to sort on

        **Returns**
            A list of intervals in the form [(strt, end), ...]
        '''
        intervals = self.chroms_list[chrom]
        # sort by the end-point
        intervals.sort(key=lambda x: (x[1], (x[1] - x[0])))   # O(n*log n)

        root = intervalnode(intervals[0][0], intervals[0][1], other=intervals[0])
        tree = reduce(lambda tree, x: tree.insert(x[0], x[1], other=x), intervals[1:], root)

        result = []

        for smallest in intervals: # O(n) (without the loop body)
            # pop the interval with the smallest end-point, keep it in the result
            if smallest.removed:
                continue # skip removed nodes
            smallest[2] = True
            result.append([smallest[0], smallest[1]]) # O(1)

            # remove (mark) intervals that overlap with the popped interval
            tree.intersect(smallest[0], smallest[1], lambda x: setattr(x.other, 'removed', True))
        return result

class intervalnode(object):
    def __init__( self, start, end, linenum=0, other=None):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority = math.ceil( (-1.0 / math.log(.5)) * math.log( -1.0 / (random.uniform(0,1) - 1)))
        self.start = start
        self.end = end
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None
        self.linenum = linenum
        self.other = other

    def insert(self, start, end, linenum=0, other=None):
        root = self
        if start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert(start, end, linenum, other)
            else:
                self.right = intervalnode(start, end, linenum, other)
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert(start, end, linenum, other)
            else:
                self.left = intervalnode(start, end, linenum, other)
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()

        if root.right and root.left:
            root.maxend = max(root.end, root.right.maxend, root.left.maxend)
            root.minend = min(root.end, root.right.minend, root.left.minend)
        elif root.right:
            root.maxend = max(root.end, root.right.maxend)
            root.minend = min(root.end, root.right.minend)
        elif root.left:
            root.maxend = max(root.end, root.left.maxend)
            root.minend = min(root.end, root.left.minend)
        return root

    def rotateright(self):
        root = self.left
        self.left = self.left.right
        root.right = self

        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)

        return root

    def rotateleft( self ):

        root = self.right
        self.right = self.right.left
        root.left = self

        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)

        return root

    def intersect( self, start, end, report_func ):
        if start < self.end and end > self.start:
            report_func(self)

        if self.left and start < self.left.maxend:
            self.left.intersect(start, end, report_func)

        if self.right and end > self.start:
            self.right.intersect(start, end, report_func)

    def traverse(self, func):
        if self.left:
            self.left.traverse(func)

        func(self)

        if self.right:
            self.right.traverse(func)

if __name__ == "__main__":
    def test_func(node):
        print("[%d, %d), %d" % (node.start, node.end, node.maxend))

    def bad_sect(lst, int_start, int_end):
        return [
            (start, end)
            for start, end in lst
            if int_start < end and int_end > start
        ]

    test = None
    intlist = []
    for x in range(20000):
        start = random.randint(0,1000000)
        end = start + random.randint(1, 1000)
        if test:
            test = test.insert(start, end)
        else:
            test = intervalnode(start, end)
        intlist.append( (start, end) )

    starttime = time.process_time()
    for x in range(50000):
        start = random.randint(0, 10000000)
        end = start + random.randint(1, 1000)
        result = []
        test.intersect(start, end, lambda x: result.append(x.linenum))
    print("%f for tree method" % (time.process_time() - starttime))

    starttime = time.process_time()
    for x in range(50000):
        start = random.randint(0, 10000000)
        end = start + random.randint(1, 1000)
        bad_sect(intlist, start, end)
    print("%f for linear (bad) method" % (time.process_time() - starttime))
