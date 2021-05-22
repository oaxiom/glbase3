
"""

Experimental-SOM like but completely pixel-to-gene clustering

"""

import random, math
import operator
import numpy
import matplotlib.pyplot as plot
from . import config

valid_init_methods = set(["random", ])
valid_metrics = set(["sum_of_deltas", ]) #

class rigidgrid:
    def __init__(self, numpy_data):
        """
        **Purpose**
            Setup the rigid grid
            
        **Arguments**
            numpy_data
                Table of expression data, rows are genes (pixels), columns are samples
        
        """
    
        self.data = numpy_data
        
    def init(self, method="random", iterations=10, metric='sum_of_deltas', seed=1234):
        """
        **Purpose**
            Initialise the starting setting for the RG.
            
        **Arguments**
            method (Optional, default=random)
                Specify the method to determine the initial starting conditions
                
                Valid methods:
                    random - randomly arrange the points in the grid.
            
            seed (Optional, default=1234)
                seed for the random number generator (if required)
            
            iterations (Optional, default=1000)
                number of iterations to plow through
                
            metric (Optional, default='quantization_error')
                metric to use to compare between two pixels for similarity.
                
                Valid metrics:
                    stdev - Sum of the differences.
                           
        """
        assert method in valid_init_methods, "'%s' method not found" % method
        assert metric in valid_metrics, "'%s' metric not found" % metric
        
        random.seed(seed)
        
        self.method = method
        self.iterations = iterations
        self.metric = metric
        
        self.sq = int(math.ceil(math.sqrt(self.data.shape[0]))) # The size of the square grid required
        self.samples_required = self.sq * self.sq
        
        config.log.info("rigidgrid(): grid set to %s*%s" % (self.sq, self.sq))
        
        self.grid = numpy.zeros([self.sq, self.sq], dtype=int) # idxs of current genes
        
        if len(self.data) < self.samples_required:
            to_add = self.samples_required - len(self.data) 
            config.log.info("rigidgrid(): I need to add %s pseudo samples to form a perfect square" % to_add)
            self.data = numpy.vstack((self.data, numpy.zeros([to_add, self.data.shape[1]])))
            
        # Init the grid
        func_init_method_map = {"random": self.__init_random}
        func_init_method_map[self.method]()
    
        config.log.info("rigidgrid(): Done init")
        
    def run(self, progress=True, verbose=True, draw_view_each_iter=False):
        """
        **Purpose**
            Run the simulation
            
        **Arguments**
            verbose (Optional, default=False)
                verbosely report what's going on.
                
            progress (Optional, default=True)
                report the current progress.
                
            draw_view_each_iter (Optional, default=False)
                draw a view of the data at each iteration
                
                (WARNING: SLOW!)
        
        **Returns**
            None
            
        """
        if draw_view_each_iter:
            self.view("debug_i0.png")

        func_metrics_map = {"sum_of_deltas": self._metric_sum_of_deltas}

        moved = 0

        # You need to pick the grid location properly...

        for iter in range(self.iterations):
            for x in range(2, self.sq-2):
                for y in range(2, self.sq-2):
                    thisi = int(x / self.sq) + int(y % self.sq)
                    thisp = self.data[thisi,:]

                    # test against all neighbours
                    near_neighbours = []
                    for xx in [-1, 0, 1]:
                        for yy in [-1, 0, 1]:
                            if xx != 0 or yy != 0:
                                thati = int((x+xx) / self.sq) + int((y+yy) % self.sq)
                                thatp = self.data[thati,:]
                                near_neighbours.append((thatp, x+xx, y+yy, func_metrics_map[self.metric](thisp, thatp))) # Yeah, sue me. 
                                print(x, y, x+xx, y+yy, thisp, thatp, near_neighbours[-1])

                    far_neighbours = []
                    for xx in [-2, 0, 2]:
                        for yy in [-2, 0, 2]:
                            if xx != 0 or yy != 0:
                                thati = int((x+xx) / self.sq) + int((y+yy) % self.sq)
                                thatp = self.data[thati,:]
                                far_neighbours.append((thatp, x+xx, y+yy, func_metrics_map[self.metric](thisp, thatp))) # Yeah, sue me. 
                                print(x, y, x+xx, y+yy, thisp, thatp, far_neighbours[-1])
                    # find the best match
                    near_neighbours = sorted(near_neighbours, key=operator.itemgetter(3))
                    far_neighbours = sorted(far_neighbours, key=operator.itemgetter(3))

                    if far_neighbours[3] > near_neighbours[3]:
                        print("Hit")
                    #print tests
                    best = tests[0]
                    # best test so swap:
                    print(best[3])
                    if best[3] > 1e-6: # minimum thrreshold to stop thrashing
                        print("Moved %s -> %s" % ("%s,%s" % (x, y), "%s,%s" % (best[1], best[2])))
                        self.grid[x,y], self.grid[best[1], best[2]] = self.grid[best[1], best[2]], self.grid[x,y] # swap
                        moved += 1
            if draw_view_each_iter:
                self.view("debug_i%s.png" % iter)

            config.log.info("rigidgrid(): '%s' iteration, moved %s items" % (iter, moved))
            moved = 0
        
    def view(self, filename, text=False, smoothed=False, **kargs):
        """
        **Purpose**
            Draw a view of the current data
        """
        dim = self.data.shape[1]

        no_row_in_plot = dim/20 + 1 #6 is arbitrarily selected
        no_col_in_plot = dim if no_row_in_plot <=1 else 20
        h = 0.2
        w = 0.001

        fig = plot.figure(figsize=(no_col_in_plot*1.5*(1+w), no_row_in_plot*1.5*(1+h)))

        for axis_num in range(self.data.shape[1]):       
            ax = fig.add_subplot(no_row_in_plot, no_col_in_plot, axis_num+1)

            this_grid = self.data[self.grid.ravel(),axis_num].reshape(self.sq, self.sq) # rearranged

            if not smoothed:        
                pl = plot.pcolormesh(this_grid[::-1], cmap='RdBu_r', vmin=self.data.min(), vmax=self.data.max())
            else:
                plot.imshow(this_grid, cmap='RdBu_r', vmin=self.data.min(), vmax=self.data.max())

            plot.axis('off')

            if text:
                print(ind, self.compname[ind])
                plt.title(self.compname[ind])
                font = {'size': text_size}

            ax.tick_params(top="off", bottom="off", left="off", right="off")
            ax.set_yticklabels("")
            ax.set_xticklabels("")
            ax.grid(False)

        plot.subplots_adjust(hspace=h*10,wspace=w*10)

        fig.savefig(filename, bbox_inches='tight', transparent=False)
        plot.close(fig)

    # ----------- Init methods

    def __init_random(self):
        """
        (Internal)
        
        initialise grid indeces randomly
        """
        # get an array of shuffled indeces:
        shuffled = list(range(0, self.samples_required)) # Not yet
        random.shuffle(shuffled) # Now it is
        for i, r in enumerate(shuffled):
            x = int(i / self.sq)
            y = int(i % self.sq)
            self.grid[x, y] = r

    # ------------ Metric methods:
    
    def _metric_sum_of_deltas(self, a, b):
        """
        (Internal)
        """
        #am = numpy.mean(a)
        #bm = numpy.mean(b)
        r = a - b
        #print a - am, b - bm, r, sum(r)
        return(abs(numpy.mean(r)))
        
if __name__ == "__main__":
    r = rigidgrid([[1,2,3],[2,3,4]])
    print(r._metric_sum_of_deltas(numpy.array([  22.,   2. ,  3.  , 2.  ,10.]), numpy.array([ 2. , 2. , 3. , 2. , 3.]) ))