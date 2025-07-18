
import copy
import pickle
import re
from shlex import split as shlexsplit

from . import utils
from . import config
from .helpers import *
from .errors import AssertionError, UnRecognisedCSVFormatError, UnrecognisedFileFormatError, ArgumentError

class _base_genelist:
    def __init__(self):
        """
        (Internal)
        This is the base derived class for all genelists.
        It contains methods available to all implementations of genelist.
        """
        self.name = None
        self.linearData = None

    def __repr__(self):
        return "<base genelist class>"

    def __in__(self, key):
        """
        (Override)

        Confer:
        if "key" in genelist:
        """
        return key in list(self.keys())

    def __bool__(self):
        """
        Fixes:
        if genelist: # contains something
            True

        and fixes:

        len(genelist) = 0
        if genelist: # Would pass even if the genelist is empty
            False

        """
        return len(self) > 0

    #def __copy__(self):
    #    raise Exception, "__copy__() is NOT supported for genelists, use gl.deepcopy() or gl.shallowcopy()"

    def __shallowcopy__(self):
        raise Exception("__shallowcopy__() is NOT supposrted for genelists, use gl.deepcopy() or gl.shallowcopy()")

    def __deepcopy__(self, fake_arg):
        raise Exception("__deepcopy__() is NOT supported for genelists, use gl.deepcopy() or gl.shallowcopy()")

    def deepcopy(self):
        """
        Confer copy to mean a deepcopy as opposed to a shallowcopy.

        This is required as genelists are compound lists.
        """
        return pickle.loads(pickle.dumps(self, -1)) # This is 2-3x faster and presumably uses less memory

    def shallowcopy(self):
        """
        (New)

        Some weird behaviour here, I know, this is so I can still get access to
        the shallow copy mechanism even though 90% of the operations are copies.
        """
        return copy.copy(self) # But doesnt this just call __copy__() anyway?

    def __len__(self):
        """
        (Override)
        get the length of the list
        """
        return len(self.linearData)

    def __int__(self):
        """
        (Override)
        get the length of the list
        NOTE: It's possible this is a bug/feature.
        I don't remove it at the moment as I'm not sure if it is used anywhere.

        """
        return len(self.linearData)

    def __iter__(self):
        """
        (Override)
        make the genelist behave like a normal iterator (list)
        """
        return self.linearData.__iter__()

    def __getitem__(self, index):
        """
        (Override)
        confers a = geneList[0] behaviour

        This is a very slow way to access the data, and may be a little inconsistent in the things
        it returns.

        NOTE:
        a = genelist[0] # returns a single dict
        a = genelist[0:10] # returns a new 10 item normal python list.
        a = genelist["name"] returns a python list containing a vertical slice of all of the "name" keys

        """
        newl = False
        if isinstance(index, int):
            # this should return a single dictionary.
            return self.linearData[index]
        elif isinstance(index, str):
            # returns all labels with that item.
            return self._findAllLabelsByKey(index)
        elif isinstance(index, slice):
            # returns a new genelist corresponding to the slice.
            newl = self.shallowcopy()
            newl.linearData = utils.qdeepcopy(self.linearData[index]) # separate the data so it can be modified.
            newl._optimiseData()
        return newl # deep copy the slice.

    def __setitem__(self, index, *args):
        """
        (Override)
        Block key editing.
        """
        raise AssertionError("Cannot modify list in-place")

    def __hash__(self):
        """
        (Override)

        compute a sensible hash value
        """
        try:
            return hash(self.name + str(self[0]) + str(self[-1]) + str(len(self))) # hash data for comparison.
        except Exception:
            try:
                return hash(self.name + str(self[0]) + str(self[-1])) # len() probably not available (delayedlist?).
            except Exception: # I bet the list is empty.
                return hash(self.name)

    def __and__(self, gene_list):
        """
        (Override)
        confer and like behaviour: c = a & b
        """
        if not self.__eq__(gene_list):
            return geneList() # returns an empty list.

        newl = self.shallowcopy()
        newl.linearData = []
        for item1 in self.linearData:
            for item2 in gene_list.linearData:
                if item1 == item2:
                    newl.linearData.append(copy.deepcopy(item1))

        newl._optimiseData()
        return newl

    def __or__(self, gene_list):
        """
        (Override)
        confer append like behaviour: c = a | b
        OR does not keep duplicates.
        """
        if not self.__eq__(gene_list):
            return geneList()
        newl = self.deepcopy()
        alist = self.linearData + gene_list.linearData
        # remove conserved duplicates;
        ulist = []
        newl = self.shallowcopy()
        for item in alist:
            if item not in ulist:
                ulist.append(item)
                newl.linearData.append(copy.deepcopy(item))

        newl._optimiseData()
        return newl

    def __add__(self, gene_list):
        """
        (Override)
        confer append like behaviour: c = a + b
        keeps duplicates (just concatenate's lists)
        """
        # check they are the same type:
        assert type(self) == type(gene_list), '"+" only works with identical types'

        # Shorcut if the genelist is empty
        if len(self.linearData) == 0:
            newl = utils.qdeepcopy(self)
            newl.linearData = utils.qdeepcopy(gene_list.linearData)
            newl._optimiseData()
            return newl

        mkeys = self._collectIdenticalKeys(gene_list)
        if not mkeys: # unable to match.
            config.log.error("No matching keys, the resulting list would be meaningless")
            return False

        newl = utils.qdeepcopy(self)
        newl.linearData.extend(utils.qdeepcopy(gene_list.linearData))
        newl._optimiseData()
        return newl

    def __sub__(self, gene_list):
        """
        (Override)
        confer c = a - b ability.
        Actually xor?
        """
        mkeys = self._collectIdenticalKeys(gene_list)
        if not mkeys: # unable to match.
            config.warning("Warning: No matching keys, unable to perform subtraction")
            return False

        newl = self.shallowcopy()
        newl.linearData = []

        dontAdd = False
        for item in self.linearData: # do a map here...
            for item2 in gene_list.linearData:
                for k in mkeys:
                    if item[k] == item2[k]:
                        dontAdd = True
                    else:
                        dontAdd = False # all mkeys must match
            if not dontAdd:
                newl.linearData.append(copy.deepcopy(item))
            dontAdd = False
        newl._optimiseData()
        return newl

    def __eq__(self, gene_list):
        """
        (Internal)
        Are the lists equivalent?
        lists now, must only have one identical key.

        This is just testing the keys...
        Wrong...
        """
        # check the hash's first to see if they are identical.
        # This is diabled as it can be very slow.
        #if self.__hash__() == gene_list.__hash__():
        #    return True

        for key in self.linearData[0]:
            if key in gene_list.linearData[0]:
                return True # just one key in common required.
        return False

    def __ne__(self, gene_list):
        """
        (Internal)
        Are the lists equivalent?
        ie do they have the same keys?
        """
        return not self.__eq__(gene_list)

    def keys(self):
        """
        return a list of all the valid keys for this geneList
        """
        if not self.linearData: return [] # Match python dict default
        return [key for key in self.linearData[0]] # Not exhaustive

    def _guessDataType(self, value):
        """
        (Internal)

        Take a guess at the most reasonable datatype to store value as.
        returns the resulting data type based on a list of logical cooercions
        (explain as I fail each cooercion).
        Used internally in _loadCSV()
        I expect this will get larger and larger with new datatypes, so it's here as
        as a separate function.

        Datatype coercion preference:
        float > list > int > location > string
        """

        try: # see if the element is a float()
            if "." in value: # if no decimal point, prefer to save as a int.
                return float(value)
            elif 'e' in value: # See if we can coocere from scientific notation
                return float(value)
            else:
                raise ValueError

        except ValueError:
            try:
                # Potential error here if it is a list of strings?
                if '[' in value and ']' in value and ',' in value and '.' in value: # Probably a Python list of floats
                    return [float(i) for i in value.strip(']').strip('[').split(',')]
                elif '[' in value and ']' in value and ',' in value: # Probably a Python list of ints
                    return [int(i) for i in value.strip(']').strip('[').split(',')]
                else:
                    raise ValueError

            except ValueError:
                try: # see if it's actually an int?
                    return int(value)
                except ValueError:
                    try: # see if I can cooerce it into a location:
                        # Turns out ~12% of loading was spent in this test:
                        if ':' in value and '-' in value:
                            return location(loc=value)
                        else:
                            raise ValueError
                    except (TypeError, IndexError, AttributeError, AssertionError, ValueError): # this is not working, just store it as a string
                        return str(value).strip()

        return "" # return an empty datatype.
        # I think it is possible to get here. If the exception at int() or float() returns something other than a
        # ValueError (Unlikely, Impossible?)

    def _processKey(self, format, column):
        """
        (Internal)
        the inner part of _loadCSV() to determine what to do with the key.
        Better in here too for security.
        """

        d = {}
        for key in format:
            if key not in ignorekeys: # ignore these tags
                #if not key in d:
                #    d[key] = {}
                if '__ignore_empty_columns' in format and format['__ignore_empty_columns']:
                    # check the column exists, if not, pad in an empty value
                    try:
                        column[format[key]]
                    except IndexError:
                        d[key] = '' # Better than None for downstream compatability
                        continue

                if isinstance(format[key], dict) and "code" in format[key]:
                    # a code block insertion goes here - any valid lib and one line python code fragment
                    # store it as a dict with the key "code"
                    d[key] = eval(format[key]["code"])
                elif isinstance(format[key], str) and "location" in format[key]:
                    # locations are very common, add support for them out of the box:
                    d[key] = eval(format[key])
                else:
                    d[key] = self._guessDataType(column[format[key]])

            elif key == "gtf_decorators": # special exceptions for gtf files
                gtf = column[format["gtf_decorators"]].strip()
                for item in gtf.split("; "):
                    if item:
                        item = item.strip()
                        ss = shlexsplit(item)
                        key = ss[0]
                        value = ss[1].strip('"')
                        d[key] = self._guessDataType(value)
        return d

    def save(self, filename=None):
        """
        **Purpose**

            Save the genelist as a binary representation.
            This is guaranteed to be available for all geneList representations, with
            the only exception being the delayedlists. As that wouldn't
            make any sense as delayedlists are not copied into memory.

            You can use this method to cache the file. It's particularly useful for large files
            that get processed once but are then used a lot.

            loading the list back into memory is relatively quick.

            list = glload("path/to/filename.glb")

            I generally used extension is glb. Although you can use
            whatever you like.

        **Arguments**

            filename
                filename (and path, if you like) to save the file to

        **Result**

            returns None
            Saves a binary representation of the geneList

        """
        assert filename, "no filename specified"

        with open(filename, "wb") as oh:
            pickle.dump(self, oh, -1)
        config.log.info(f"Saved binary version: '{filename}'")

    def from_pandas(self, pandas_data_frame):
        """
        **Purpose**

            Convert a pandas dataFrame to a genelist

            NOTE: This is an INPLACE method that will REPLACE any exisiting data
            in the genelist

        **Arguments**

            pandas_data_frame (Required)
                The pandas data frame to convert

        **Result**
            None
            The object is populated by

        """
        if len(self) > 0:
            config.log.warning('genelist.from_pandas() will overwrite the existing data in the genelist')

        newl = []
        key_names = pandas_data_frame.columns
        for index, row in pandas_data_frame.iterrows():
            newitem = {k: item for k, item in zip(key_names, row)}
            newl.append(newitem)
        self.linearData = newl
        self._optimiseData()

        config.log.info("genelist.from_pandas() imported dataFrame")
