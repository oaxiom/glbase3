"""

Init wrapper around PEBL 1.0.2.

I will gradually port functionality into this module as time goes by.

Part of glbase, but not copyright by me.

"""

import tempfile

import os, sys, numpy
from bisect import insort, bisect
import networkx as nx
try:
	from pygraphviz import AGraph
except Exception:
	try:
		from graphviz import AGraph # Throw an error if failed here
	except Exception:
		pass # This is pretty bad I think...

import matplotlib.pyplot as plot

from . import data, evaluator, network, posterior, prior, result 
from . import learner.greedy, learner.custom, learner.simanneal, learner.exhaustive 


from .. import config
from .. import progress

valid_learners = set(["greedy", "exhaustive"])

# helpers for consensus_network() move to bayesutils?
def colorize_edge(weight):
    colors = "dcba986420"
    breakpoints = [.1, .2, .3, .4, .5, .6, .7, .8, .9]
    return "#" + str(colors[bisect(breakpoints, weight)])*6

def node(n, position):
    s = "\t\"%s\"" % n.name
    if position:
        x,y = position
        s += " [pos=\"%d,%d\"]" % (x,y)
    return s + ";"

class bayes:
    def __init__(self, parent):
        """
        **Purpose**
            Wrapper around bayes
        
            A typical entry should not be direct, it should be something like:
            
            expn.bayes.learn()
            expn.bayes.run()
            ...
        
        **Arguments**
            parent
                The parent genelist
                
        """        
        self.parent = parent
        self.learner = None
        self.__dataset = data.Dataset(self.parent.numpy_array_all_data.T, 
            samples=numpy.array(self.parent.getConditionNames()), 
            variables=numpy.array([data.Variable(i) for i in numpy.array(self.parent["name"])])) # Just a hack for now
        self.__dataset.discretize(numbins=10) # required, although numbins is optional
    
    def learn(self, mode="greedy", num_learners=1, max_iterations=10000):
        """
        **Purpose**
            Do Bayes learning
        
            By default learn is set to run a quick greedy learner for testing purposes.
            You will need to play with the settings for a real run.
        
        **Arguments**
            mode (Optional, default="greedy")
                the learner method to use, 
                valid learners are %s
                
            num_learners (Optional, default=1)
                The number of learners to use
            
            max_iterations (Optional, default=10000)
                number of iterations to use in each learner.
                
                
        """ % valid_learners
        assert self.__dataset, "bayes: init failed..."
        assert mode in valid_learners, "'%s' learner not found" % mode
        
        if mode == "greedy":
            self.learners = [learner.greedy.GreedyLearner(self.__dataset, max_iterations=max_iterations) for i in range(num_learners)]
        elif mode == "exhaustive": 
            self.learners = [learner.greedy.ListLearner(self.__dataset) for i in range(num_learners)]
            
        config.log.info("bayes.learn: Set %s %s learners" % (num_learners, mode))

    def run(self):
        assert self.learners, "bayes: No learner, do learn() first"
        
        runs = []
        p = progress.progressbar(len(self.learners))
        for i, l in enumerate(self.learners):
            runs.append(l.run())
            p.update(i)
        self.result = result.merge(runs)
        config.log.info("bayes.run: Finished")
    
    def results(self, dirname):
        self.result.tohtml(dirname)
        config.log.info("bayes.results: Wrote '%s'" % dirname) 
        
    def consensus_network(self, filename, threshold, write_gml=False):
        """
        **Purpose**
            Draw a consensus network for a specified threshold
            
        **Arguments**
            filename (Required)
                filename to save the png and svg to. 
                
            threshold (Required)
                threshold value 0..1.0 for the consensus network
                
            write_gml (Optional, default=False)
                write the GML file for loadinging into e.g. Cytoscape
        
        
        """
        assert self.result, "bayes.consensus_network: Bayes network not learnt yet"
        assert filename, "You must specify a filename"
        assert threshold, "You must specify a threshold"
        assert threshold <= 1.0, "threshold is >= 1.0!"
        assert threshold >= 0.0, "threshold is <= 0.0!"
                   
        #print [n.score for n in list(reversed(self.result.networks))]

        this_posterior = posterior.from_sorted_scored_networks(self.result.nodes, list(reversed(self.result.networks)))
        net_this_threshold = this_posterior.consensus_network(threshold)

        cm = this_posterior.consensus_matrix # A method/property...
        top = this_posterior[0] # posterior object
        top.layout()
        nodes = net_this_threshold.nodes

        A = AGraph()
        net = A.from_string(string=net_this_threshold.as_dotstring()) # Load their Network into a AGraph and NX.
        G = nx.from_agraph(net)
        edge_cols = [colorize_edge(cm[src][dest]) for src, dest in net_this_threshold.edges]
        #pos = nx.graphviz_layout(G, 'dot')
        #fig = plot.figure()
        #ax = fig.add_subplot(111)
        #nx.draw_graphviz(G, 'dot', with_labels=True, edge_color=edge_cols, linewidths=0)
        #fig.tight_layout()
        #fig.savefig(filename.replace(".png", ".pgv.png"))

        dotstr = "\n".join(
            ["digraph G {"] + 
            [node(n, pos) for n, pos in zip(nodes, top.node_positions)] + 
            ["\t\"%s\" -> \"%s\" [color=\"%s\"];" % (nodes[src].name, nodes[dest].name, colorize_edge(cm[src][dest])) for src,dest in net_this_threshold.edges] + 
            ["}"])

        fd, fname = tempfile.mkstemp()
        open(fname, 'w').write(dotstr)

        if write_gml:
            nx.write_gml(G, filename.replace(".png", ".gml"))
        
        # Better done through pygraphviz gives different results...
        os.system("dot -n1 -Tpng -o%s %s" % (filename, fname))
        os.system("dot -n1 -Tsvg -o%s %s" % (filename.replace(".png", ".svg"), fname))
        os.remove(fname)
        config.log.info("bayes.consensus_network: Saved '%s'" % filename)