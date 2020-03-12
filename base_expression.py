"""

Expression class

Basic handling for microarray and rna-seq and realtime PCR like data

"""

import sys, os, csv, string, math, collections

from operator import itemgetter

import numpy
from numpy import array, arange, meshgrid, zeros, linspace, mean, object_, std # This use of array here is not good.

from . import config
from .flags import *
from .draw import draw
from .genelist import genelist
from .progress import progressbar
from .errors import AssertionError, ArgumentError, ExpressionNonUniqueConditionNameError
from .utils import qdeepcopy

class base_expression(genelist):
    def __init__(self, filename=None, loadable_list=None, format=None, expn=None, silent:bool=False, **kargs):
        """
        See the documentation in the expression class.

        This is the underlying base expression object and is not designed for direct usage.
        """
        '''
        if not loadable_list:
            # these are only required if not loading a list
            assert expn, "'expn' argument cannot be empty"
            assert filename, "no filename to load"
            assert format, "required argument 'format' is missing"
            assert os.path.exists(os.path.realpath(filename)), "'%s' not found" % filename
        else:
            # probably should put some more sanity checking in here.
            assert loadable_list[0], "the list to load does not appear to be a proper list"
        '''

        if "cv_err" in kargs or "err_up" in kargs or "err_dn" in kargs:
            raise NotImplementedError("Whoops! I haven't finished expression class - cv_err, err_up and err_dn are not implemented")

        valig_args = ["cond_names", "name", "force_tsv", "nan_value"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError(self.__init__, k)

        genelist.__init__(self)

        self.filename = filename
        self._conditions = [] # Provide a dummy conditions temporarily
        self.name = "None"

        if "name" in kargs and kargs["name"]:
            self.name = kargs["name"]
        elif filename:
            self.name = "".join(self.filename.split(".")[:-1])

        if not loadable_list and not expn:
            config.log.info("expression: made an empty expression object")
            return()

        if loadable_list:
            self.load_list(loadable_list, expn, **kargs)
        else:
            # This is a placeholder at the moment,
            # I reload the expn and err values back into the format
            # When you redo this, remember to also redo load_list()
            newf = format
            newf["conditions"] = {"code": expn}
            if "err" in kargs and kargs["err"]:
                newf["err"] = {"code": kargs["err"]}
            elif "cv_err" in kargs and kargs["cv_err"]:
                newf["cv_err"] = kargs["cv_err"]

            if "force_tsv" in kargs and kargs["force_tsv"]:
                newf["force_tsv"] = True

            format = newf

            self.loadCSV(filename=filename, format=format) # no need for error checking here - it's in genelist now.

            if "cond_names" in kargs and kargs["cond_names"]:
                self._conditions = kargs["cond_names"]
            else:
                # re-open the file and try to guess the conditions
                # reopen the file to get the condition headers.
                oh = open(filename, "rU")
                if "force_tsv" in format and format["force_tsv"]:
                    reader = csv.reader(oh, dialect=csv.excel_tab)
                elif "dialect" in format:
                    reader = csv.reader(oh, dialect=format["dialect"])
                else:
                    reader = csv.reader(oh)

                do = False
                self._conditions = []
                for index, column in enumerate(reader):
                    if "skiptill" in kargs:
                        if kargs["skiptill"] in column:
                            do = True
                    elif "skiplines" in kargs:
                        if index == kargs["skiplines"]:
                            do = True
                    else:
                        do = True # do anyway

                    if do:
                        names = eval("{0}".format(format["conditions"]["code"])) # yay, more nice happy arbitrary code execution.

                        if names:
                            self._conditions = [str(k) for k in names]
                        break
                oh.close()

                if not silent:
                    config.log.info("expression: I found the following conditions:")
                    config.log.info("\n".join(["%s\t%s" % (n, i) for n, i in enumerate(self._conditions)]))

        # coerce the conditions errs etc to floats
        nans = set(('nan', 'Nan', 'NaN'))
        for idx, i in enumerate(self):
            try:
                # Nan policy:
                if True in [t in nans for t in i["conditions"]]:
                    config.log.warning("line {0}, contains Nan, filling with 0".format(idx))
                    newc = []
                    for c in i['conditions']:
                        if c in nans:
                            newc.append(0.0) # nan policy
                        else:
                            newc.append(c)
                    i['conditions'] = newc
                i["conditions"] = [float(str(t).replace(",", "")) for t in i["conditions"]] # because somebody once sent me a file with ',' for thousands!
            except ValueError:
                config.log.warning("line %s, contains missing data (%s), filling with 0" % (idx, i["conditions"]))
                i["conditions"] = [0 for t in self._conditions] # Use conditions as the example I had here was also missing all of the other values.

            # These will bomb on missing data...
            if "err" in i:
                i["err"] = [float(t) for t in i["err"]]
            if "cv_err" in i:
                i["cv_err"] = [float(t) for t in i["cv_err"]]

        self.__check_condition_names_are_unique()
        self._optimiseData()
        if not silent:
            config.log.info("expression: loaded %s items, %s conditions" % (len(self), len(self.getConditionNames())))

    def __check_condition_names_are_unique(self):
        """
        Bit of gotcha this one, but expression objects must have unique condition names
        or lots of things break. Here, check the condition names are unique.

        """
        if len(self._conditions) > len(set(self._conditions)):
            raise ExpressionNonUniqueConditionNameError(self._conditions)
        return(False)

    def __repr__(self):
        return("glbase.expression")

    def _load_numpy_back_into_linearData(self):
        """
        For routines that make a change in self.numpy_array_all_data

        this must be called after to propogate the changes back into linearData

        """
        for i, row in enumerate(self.numpy_array_all_data):
            self.linearData[i]["conditions"] = list(row)
        self._optimiseData()

    def _optimiseData(self):
        """
        (Override)
        (Internal)
        Add expression optimisations
        """
        genelist._optimiseData(self) # do the parent optimise.

        # generate a serialised version of the array conditions.
        self.numpy_array_all_data = numpy.array([i["conditions"] for i in self.linearData])

        # could be done with dict comp:
        data = {}
        for index, name in enumerate(self._conditions):
            if not name in data:
                data[name] = self.numpy_array_all_data[:,index]
        self.serialisedArrayDataDict = data

        # list;
        self.serialisedArrayDataList = [self.serialisedArrayDataDict[key] for key in self._conditions]
        #self.serialisedArrayDataList = all_array_data # This consumes massive amounts of memory.
        # presumably something downstream is doing something nasty.

        return(True)

    def saveCSV(self, filename=None, interleave_errors=True, no_header=False, no_col1_header=False, **kargs):
        """
        A CSV version of saveTSV(), see saveTSV() for syntax
        """
        self.saveTSV(filename=filename, tsv=False, interleave_errors=True, no_header=False, no_col1_header=False, **kargs)
        config.log.info("saveCSV(): Saved '%s'" % filename)

    def saveTSV(self, filename=None, tsv=True, interleave_errors=True, no_header=False, no_col1_header=False, **kargs):
        """
        (Override)
        **Purpose**
            Save the microarray data as a tsv file
            This is a little different from the normal genelist.saveTSV()
            as I want to make certain that the condition data is written in a sensible manner at
            the end of the TSV.
            I also need to deal with grid like structures etc.

            As a general warning, use expression.save() in preference to this.
            This save is not guaranteed to survive reloading into glbase, and is particularly
            troublesome in the case of expression objects. Indeed, the default guesser when loading
            a genelist object will incorrectly load an expression object with error values
            and will probably bodge any other arrangement too.

        **Arguments**
            filename
                The filename (with a valid path) to save the file to.

            interleave_errors (Optional, default=True)
                By default the errors are interleaved so that the sample data will be arranged:

                Sample1 Err1 Sample2 Err2

                if interleave_errors=False then:

                Sample1 Sample2 Err1 Err2

            no_col1_header (Optional, default=False)
                In case you want a table like this:

                    A   B   C   D
                W   1   2   3   4
                X   2   2   2   2
                Y   2   2   2   2
                Z   2   2   2   2

                i.e. the top left column label is empty.

        **Returns**
            returns None

        """
        self._save_TSV_CSV(filename=filename, tsv=tsv, interleave_errors=True, no_header=False, no_col1_header=False, **kargs)
        config.log.info("saveTSV(): Saved '%s'" % filename)

    def _save_TSV_CSV(self, filename=None, tsv=True, interleave_errors=True, no_header=False, no_col1_header=False, **kargs):
        """
        Internal unified saveCSV/TSV for expression objects
        """
        valig_args = ["filename", "tsv", "key_order", "no_header"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError(self.saveCSV, k)

        assert filename, "you must specify a filename"

        oh = open(os.path.realpath(filename), "w")
        if tsv:
            writer = csv.writer(oh, dialect=csv.excel_tab)
        else:
            writer = csv.writer(oh)

        array_data_keys = ("conditions", "err", "cv_err")

        write_keys = []
        if "key_order" in kargs:
            write_keys = kargs["key_order"]
            # now add in any missing keys to the right side of the list:
            for item in list(self.keys()):
                if item not in write_keys and item not in array_data_keys: # But omit the array_data_keys
                    write_keys.append(item)
        else:
            # just select them all:
            write_keys = [k for k in list(self.keys()) if not k in array_data_keys]

        if "err" in list(self.keys()):
            if interleave_errors:
                conds = ["mean_%s" % c for c in self.getConditionNames()]
                errs = ["err_%s" % c for c in self.getConditionNames()]
                paired = [val for pair in zip(conds, errs) for val in pair]

                if not no_header:
                    title_row = [k for k in write_keys if k in list(self.keys())]
                    writer.writerow(title_row + paired)

                for data in self.linearData:
                    line = [data[k] for k in write_keys if k in data]

                    interleaved_data = [val for pair in zip(data["conditions"], data["err"]) for val in pair] # I never understand how these work, but what the hell.

                    writer.writerow(line + interleaved_data)# conditions go last.
                oh.close()
            else:
                if not no_header:
                    title_row = [k for k in write_keys in k in list(self.keys())]
                    writer.writerow(write_keys + self.getConditionNames() + ["err_%s" % i for i in self.getConditionNames()])

                for data in self.linearData:
                    line = [data[k] for k in write_keys if k in data]
                    writer.writerow(line + data["conditions"] + data["err"])# conditions go last.
                oh.close()

        else: # no error, very easy:
            if not no_header:
                title_row = [k for k in write_keys if k in list(self.keys())]
                if no_col1_header:
                    title_row[0] = ""
                writer.writerow(title_row + self.getConditionNames())

            for data in self.linearData:
                line = [data[k] for k in write_keys if k in data]
                writer.writerow(line + data["conditions"])# conditions go last.
            oh.close()

        return(None)

    def sort(self, key, reverse=False):
        """
        This is slightly different from the vanilla genelist's sort - you can pass it the name of
        a condition. Take care to make sure the condition name is not also a valid list key.
        The algorithm searches the genelist before searching the array for your particular condition.

        Also take care with this one: It is one of the few in-place list
        modifiers.

        **Arguments**

        key
            must be a valid key in the genelist or the name of an array condition.

        reverse (Optional, default=False)
            By default the list is sorted smallest to largest.
            reverse = True sorts largest to smallest.

        **Result**

        returns True if succesful.

        returns False if not valid.
        """
        assert (key in self.linearData[0]) or key in self._conditions, "'%s' search key not found in list or array data" % key

        if key in self.linearData[0]:
            return(genelist.sort(self, key, reverse=reverse)) # use the parents sort.
        else:
            if key in self._conditions:
                name_index = self._conditions.index(key)
                self.linearData = sorted(self.linearData, key=lambda x: x["conditions"][name_index]) # the original sort() was overridden.
                if reverse:
                    self.linearData.reverse()
                self._optimiseData()
                return(True)
        return(False)

    def load_list(self, list_to_load, expn=None, name=False, cond_names=None, nan_value=0):
        """
        **Purpose**
            You've generated your own [{ ... }, { ...}] like list
            (A list of dicts) and you want to either reload it into
            a genelist-like object or load it into an empty genelist.
            This is the method to do that officially.

            This method should be used with care. Some sanity
            checking is done. But not very much.

            This load_list is modified for expression-like genelists.
            (eg. expression()). Here you can load keys into conditions based on
            their key names.

        **Arguments**
            list_to_load
                must be a list of dicts.

            expn (optional)
                A list of key names to construct the expression data from
                If not specified then it assumes your list already has a correctly formatted
                "conditions" key.

        **Returns**
            None. This is one of the few IN PLACE methods. and returns
            None.
        """
        assert list_to_load[0], "list_to_load does not appear to be a valid list"

        __nan_warnings = False
        nans = frozenset(["Inf", "-Inf", "NA", "Nan", "NaN"])

        if expn:
            assert isinstance(expn, list), "'expn' must be a list of keys"
            # Bodge in a new "conditions" key:
            newl = []
            for i in list_to_load:
                new = i.copy()
                nl = [i[k] for k in expn]

                # test for Inf, -Inf, NA, NaN, etc.
                if True in [ti in nans for ti in nl]: # woah! Nan here.
                    t = []
                    for item in nl:
                        if item in nans:
                            t.append(nan_value)
                        else:
                            t.append(item)
                    nl = t
                    if not __nan_warnings:
                        __nan_warnings = True
                        config.log.warning("Expression list contains 'not a number' values, setting them to <nan_value=%s>" % nan_value)

                new["conditions"] = nl
                for k in expn:
                    del new[k]
                newl.append(new)
            self._conditions = expn
        else:
            newl = list_to_load
            if cond_names: # user sent the conditions names. Hope they are in the same order
                assert len(cond_names) == len(newl[0]["conditions"]), "cond_names is not the same length as the number of conditions"
                self._conditions = cond_names
            else:
                # conditions can get lost in a loadable list. fill in a dummy one
                if len(self._conditions) != len(newl[0]["conditions"]):
                    self._conditions = ["cond_%s" % i for i in range(len(newl[0]["conditions"]))]

        # Now call parent with new list
        genelist.load_list(self, newl, name)

    def from_pandas(self, pandas_data_frame, condition_names=None):
        """
        **Purpose**

            Convert a pandas dataFrame to a genelist

            NOTE: This is an INPLACE method that will REPLACE any exisiting data
            in the

        **Arguments**

            pandas_data_frame (Required)
                The pandas data frame to convert

            condition_names (Required)
                A list of Column names from the Pandas frame to use as expression data

        **Result**
            None
            The object is populated by the Pandas object

        """
        assert condition_names, 'You must specify condition_names'
        assert isinstance(condition_names, list), 'condition_names must be a list of colun names'
        if len(self) >0:
            config.log.warning('expression.from_pandas() will overwrite the existing data in the expression')

        newl = []
        key_names = pandas_data_frame.columns
        for index, row in pandas_data_frame.iterrows():
            newitem = {}

            # load normal keys:
            for k, item in zip(key_names, row):
                if k not in condition_names:
                    newitem[k] = item
            # load conditions, in-order:
            dict_items = dict(zip(key_names, row))
            newitem['conditions'] = [dict_items[z] for z in condition_names]

            newl.append(newitem)

        self._conditions = condition_names
        self.linearData = newl
        self._optimiseData()

        config.log.info("expression.from_pandas() imported dataFrame")

    def getConditionNames(self):
        """
        returns a list of the condition headers
        """
        return(list(self._conditions))

    def setConditionNames(self, new_cond_names):
        """
        rename the conditions names for the expression data

        THIS IS AN IN-PLACE method and returns None
        """
        assert len(new_cond_names) == len(self._conditions), "setConditionNames(): new and old condition names are different lengths (%s vs. %s)" % (len(new_cond_names), len(self._conditions))

        self.__check_condition_names_are_unique()
        self._conditions = list(new_cond_names)
        self._optimiseData()
        return(self._conditions)
