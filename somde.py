'''

Differential expression using SOMs.

'''

import numpy, copy, math
from . import config
from .som import SOM, normalize

class somde(SOM):
    def __init__(self, parent, name):
        SOM.__init__(self, parent, name)
    
    def preparse(self, Z, filename=None, threshold_value=False, digitize=False,):
        '''
        **Purpose**
            Preparse the expression data to firstly reduce the search space, and secondly to focus
            the SOMs around important genes.
        
        '''
        array_std = numpy.std(self.data_raw) # I have to go back ot the raw_data
        
        new_raw = []
        new_data = []
        newdlab = []
        for i, row in enumerate(self.data_raw):
            #print row
            m = numpy.mean(row)
            s = numpy.std(row)
            CV = s / m
            print(self.dlabel[i], CV)
            if CV > 1.0 and CV < 7.0:
                zzs = (row - m) / array_std
                print(self.dlabel[i], CV, zzs)
                new_raw.append(row)
                new_data.append(zzs)
                newdlab.append(self.dlabel[i])
        
        #self.data_raw = 
        self.dlabel = newdlab
        
        #
        # No need to normalise the zzs:
        self.data_raw = numpy.array(new_raw)
        self.data = numpy.copy(self.data_raw)
        self.data = numpy.copy(self.data_raw)
        if threshold_value:
            assert digitize, "If using threshold_values you must also use digitize"
            self.data -= threshold_value
            self.data[self.data<0.0] = 0

            if digitize == 1:
                self.data[self.data>0.0] = 1
            elif digitize:
                self.data /= self.data.max() 
                self.data *= digitize
                self.data = numpy.ceil(self.data) # always move up as the zeros have already been set
            self.data = normalize(self.data, method='var') # Don't use if using zzs
        else:
            self.data = normalize(numpy.copy(self.data_raw), method='var') # Don't use zzs
        
        #self.data = numpy.array(new_data)

        print(self.data.shape)
        print(len(self.dlabel))
    
        sq = math.ceil(math.sqrt(len(self.dlabel)))
        self.mapsize = (int(sq),int(sq))# I want to re-guess a reasonable mapsize after culling genes.
    
        # Have to repeat the last part of config()
        self.finalise_config()

        config.log.info("somde.preparse: %s rows survived preparse" % len(self.dlabel))