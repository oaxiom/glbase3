

import timeit, sys

sys.path.append('../../')

# Test for deepcopy or

print('location(loc="chr1:10000-10200")')
print(timeit.timeit('location(loc="chr1:10000-10200")', number=100000, setup='from location import location'))
print(timeit.timeit('location(loc="chr1:10000-10200")', number=100000,  setup='from location2 import location'))

print('l = l.expand(2000)')
print(timeit.timeit('_ = l1.expand(2000)', number=100000, setup='from location import location; l1 = location("chr1:10000-10200")'))
print(timeit.timeit('_ = l2.expand(2000)', number=100000, setup='from location2 import location; l2 = location("chr1:10000-10200")'))

print('l = l.pointify()')
print(timeit.timeit('_ = l1.pointify()', number=100000, setup='from location import location; l1 = location("chr1:10000-10200")'))
print(timeit.timeit('_ = l2.pointify()', number=100000, setup='from location2 import location; l2 = location("chr1:10000-10200")'))
