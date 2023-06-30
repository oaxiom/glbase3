"""

location.py

part of glbase.

This class is an internal class that implements a more convenient way to manipulate
genomic coordiantes.

TODO:
. This is really due for revision...


"""

import copy, pickle

class location:
    def __init__(self, loc=None, chr=None, left=None, right=None):
        if loc:
            if isinstance(loc, location):
                # It's actually already a loc.
                # I want to copy it and leave.
                self.chrom = loc.chrom
                self.left = loc.left
                self.right = loc.right

            elif isinstance(loc, str):
                s = loc.lower().replace(",", "") # ucsc includes commas, remove them so you can cut and paste
                t = s.split(":")
                self.chrom = t[0].strip("chr").rstrip().upper()
                self.left = int(t[1].split("-")[0])
                self.right = int(t[1].split("-")[1])
        else:
            self.chrom = str(chr).strip("chr").rstrip().upper()
            self.left  = int(left)
            self.right = int(right)

    def __eq__(self, other):
        if other:
            if isinstance(other, str):
                return (str(self) == str(other.replace(",", ""))) # use string comparison.

            # use a faster ? dict comparison, or throw an exception, as this item probably not a <location>
            if (
                self.chrom == other.chrom
                and self.left == other.left
                and self.right == other.right
            ):
                return True

        return False

    def __lt__(self, other): # Required for sorting
        # Make locations sortable
        if self.chrom < other.chrom:
            return True

        elif self.chrom == other.chrom:
            if self.left < other.left:
                return True
            return False
        return False

    def __hash__(self):
        return hash(self.__str__())

    def __deepcopy__(self, memo):
        return pickle.loads(pickle.dumps(self, -1)) # This is 2-3x faster and presumably uses less memory

    def __bool__(self):
        return True

    def __repr__(self):
        return f"chr{self.chrom}:{self.left}-{self.right}"

    def __len__(self):
        # work out the span.
        return max([0, self.right - self.left])

    def split(self, value=None):
        # ignores the 'value' argument completely and returns a three-ple
        return (self.chrom, self.left, self.right)

    def __getitem__(self, key):
        # Emulate a dict;
        if key == "chr":
            return self.chrom
        elif key == "left":
            return self.left
        elif key == 'right':
            return self.right
        return None

    def __setitem__(self, key, value):
        self.loc[key] = value
        self.__update()

    def __str__(self):
        return f"chr{self.chrom}:{self.left}-{self.right}"

    """
    these methods below should copy the location and send a modified version back.
    """
    def expand(self, base_pairs):
        return location(chr=self.chrom, left=self.left - base_pairs, right=self.right + base_pairs)

    def expandLeft(self, base_pairs):
        return location(chr=self.chrom, left=self.left - base_pairs, right=self.right)

    def expandRight(self, base_pairs):
        return location(chr=self.chrom, left=self.left, right=self.right + base_pairs)

    def shrink(self, base_pairs):
        return location(chr=self.chrom, left=self.left + base_pairs, right=self.right - base_pairs)

    def shrinkLeft(self, base_pairs):
        return location(chr=self.chrom, left=self.left - base_pairs, right=self.right )

    def shrinkRight(self, base_pairs):
        return location(chr=self.chrom, left=self.left, right=self.right - base_pairs)

    def pointLeft(self):
        return location(chr=self.chrom, left=self.left, right=self.left)

    def pointRight(self):
        return location(chr=self.chrom, left=self.right, right=self.right)

    def pointify(self):
        centre = (self.left + self.right) // 2
        return location(chr=self.chrom, left=centre, right=centre)

    def collide(self, loc):
        if loc.chrom != self.chrom: return False
        return self.right >= loc.left and self.left <= loc.right

    def qcollide(self, loc):
        """
        **Purpose**
            perform a collision with another location object.
            This assumes you have already checked the locations are on the same chromosome.

        **Returns**
            True or False
        """
        return self.right >= loc.left and self.left <= loc.right # nice one-liner

    def distance(self, loc):
        """
        **Purpose**
            calculate the distance between two locations.

        **Returns**
            an integer indicating the distance, note that
            the chromosomes should be the same or it will raise an
            exception. distance() should not be used as a test for
            overlap. use collide() for that.
        """
        assert self.chrom == loc.chrom, "chromosomes are not the same, {self} vs {loc}"
        return self.qdistance(loc)

    def qdistance(self, loc):
        """
        (Internal)
        ignore the assert.
        """
        return ((self.left + self.right) // 2) - ((loc.left + loc.right) // 2)

    def __sub__(self, loc):
        """
        **Purpose**
            Allow things like:

            distance = locA - locB
        """
        return self.distance(loc)

    '''
    def offset(self, base_pairs):
        """
        get a new location offset from the 5' end by n base pairs
        returns a point location.
        """
        new = copy.deepcopy(self)
        new.loc["left"] += base_pairs
        new.loc["right"] = new.loc["left"]
        new.__update()
        return location(chr=self.chrom, left=self.left + base_pairs right=self.left)
    '''

    def keys(self):
        return ('chr', 'left', 'right')

    def values(self):
        return (self.chrom, self.left, self.right)

if __name__ == "__main__":
    import timeit

    s = "a = location(loc='chr1:1000-2000').pointify()"
    t = timeit.Timer(s, "from location import location")
    print(s)
    print("%.2f usec/pass" % (1000000 * t.timeit(number=100000)/100000))

    s = "a = location(chr='chr1', left=1000, right=2000).pointify()"
    t = timeit.Timer(s, "from location import location")
    print(s)
    print("%.2f usec/pass" % (1000000 * t.timeit(number=100000)/100000))

    s = "a = location(chr='chr1', left=1000, right=2000).expand(1000)"
    t = timeit.Timer(s, "from location import location")
    print(s)
    print("%.2f usec/pass" % (1000000 * t.timeit(number=100000)/100000))

    s = "a = str(location(chr='chr1', left=1000, right=2000))"
    t = timeit.Timer(s, "from location import location")
    print(s)
    print("%.2f usec/pass" % (1000000 * t.timeit(number=100000)/100000))

    s = "a.qcollide(b)"
    t = timeit.Timer(s, "from location import location ; a = location(chr='chr1', left=1000, right=2000); b = location(chr='chr1', left=1200, right=2200)")
    print(s)
    print("%.2f usec/pass" % (1000000 * t.timeit(number=100000)/100000))
