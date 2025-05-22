"""

Multi-dimensional PCA squish



"""

import operator
import heapq
import itertools
import math
import collections
import matplotlib
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from operator import itemgetter

import numpy
from sklearn.decomposition import PCA
import networkx as nx # Should only be possible to get here if networkx is available.
import scipy.cluster.vq
import scipy.stats
import scipy.cluster.hierarchy
import scipy.spatial.distance

from . import network_support
from . import utils, config
from .progress import progressbar
from functools import reduce


class manifold_mdsquish:
    def __init__(self, parent):
        """
        **Purpose**
            Dimensionally squish complex data into 2D space, draw a tree and network, etc...
            
            You don't need to call me directly, use expn.mdsquish.squish()
        
        """       
        self.parent = parent
        self.names = parent.getConditionNames()
        self.data = None
        self.G = None # Network
        self.path_func_mapper = network_support.path_func_mapper # For checking valid paths
        self.__layout_data = None
        self.__last_thresholds = (-1,-1,-1)
        
    def squish(self, pcas=None, whiten=False, traversal_weight=1.0, mean_subtraction=False, 
        normalise_scales=True, invert_scales=True, log_transform=False, scale_by_percent_variance=False):
        """
        **Purpose**
            Calculate the mdsquish, using the specified PCAs
            
        **Arguments**
            pcas (Required)
                The range of PCAs to use for the squish.
                
                You should probably explore the most useful PCAs before deciding on the PCAs to use,
                mdsquish uses the same PCA system as expression.pca (it does not use svd)
                
                PCAs begin at 1, and not at 0 as might be expected.
                
            log_transform (Optinal, default=False)
                log transform before normalising between 0 and 1. This has the effect of flattening the curve out.
                Check the output from hist() to make sure this is doing what you expect.
            
            scale_by_percent_variance (Optional, default=False)
                scale each PC so that it is scaled to its percent variance. This leads to 
                reduced influence on the mdsquish plot from lower PCs.
                
                It's not clear exactly why, but you are advised to leave this as False
            
            normalise_scales (Optional, default=True)
                normalise the euclidean distances to 0 ... 1.0 range
                
            invert_scales (Optional, default=True)
                normally the lower the score the closer the relationship, if this is true the
                scales are inverted to the more natural:
                
                (distant) 0.0 ... ... 1.0 (close)

            whiten (Optional, default=False)
                'Normalise' the data. Each feature is divided by its standard deviation 
                across all observations to give it unit variance.
                If your samples have a lot of really big outliers this may help.
    
            traversal_weight (Optional, default=1.0)
                This specifies the cost/weight of moving from one cell to another in the network.
                Use it as a penalty for excessive cell-cell movement when you are path finding
    
        """
        assert pcas, "You must specify the PCA dimensions to use for the squish"
        
        self.matrix = self.parent.getExpressionTable().T
        
        self.__model = PCA(n_components=max(pcas)+1, whiten=whiten)
        self.__model.fit(self.matrix)
        self.__transform = self.__model.transform(self.matrix) # project the data into the PCA
        self.__explained_variance = numpy.array(self.__model.explained_variance_ratio_) * 100.0
        self.data = True
        
        # I still need to collect the appropriate PCAs incase user skips some PCs:  
                     
        pcs = []
        for pc in pcas:
            if scale_by_percent_variance:
                pcs.append(self.__transform[:,pc-1] * self.__explained_variance[pc-1])
            else:
                pcs.append(self.__transform[:,pc-1])
        pcs = numpy.array(pcs)
        
        #print pcs
        
        res = numpy.zeros((pcs.shape[1], pcs.shape[1]))
        
        # Iterate over dimensions
        # TODO : Scale the data by percent variance.
        for item1 in range(pcs.shape[1]): # each row is a PC, this is iterating the columns (samples)
            for item2 in range(pcs.shape[1]):
                if item1 != item2:
                    sample1 = pcs[:,item1]
                    sample2 = pcs[:,item2]
                    #print "s", sample1, sample2, [(i[0], i[1]) for i in zip(sample1, sample2)]
                    deltas = [abs(i[0] - i[1])**2 for i in zip(sample1, sample2)] # This is sqeuclidean...
                    euclidean_distance = reduce(operator.add, deltas, 0)
                    #print self.names[item1], self.names[item2], 1-euclidean_distance, [((i[0], i[1])) for i in (sample1, sample2)]
                    res[item1, item2] = math.sqrt(euclidean_distance) # now is actual euclidean distance
        
        if log_transform:
            res = numpy.log10(res+0.1)
        
        if normalise_scales: # Shouldn't this be column-wise?
            res = res - res.min()
            res = res / res.max()

        if invert_scales:
            res = 1 - res

        self.traversal_weight = traversal_weight

        self.distances = res
        
        config.log.info("mdsquish.squish: done squishing %s..." % pcas)
    
    def hist(self, filename, **kargs):
        """
        **Purpose**
            save a histogram of the distances (KDE smoothed)
        
        **Arguments**
            filename (Required)
                the filename to save to.
        """
        fig = self.parent.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        # I just want to flatten all distances
        d = self.distances.flatten()
        values = utils.kde(d, range=(self.distances.min(), self.distances.max()), covariance=0.1, bins=1000)
        ax.plot(values, color="grey")
        
        ax.set_xticks(list(range(0, 1000, 50)))
        ax.set_xticklabels(numpy.arange(self.distances.min(), self.distances.max(), self.distances.max() * (50.0/1000.0)))
        ax.grid(True)
    
        self.parent.draw.do_common_args(ax, **kargs)
        real_filename = self.parent.draw.savefigure(fig, filename)
        
        config.log.info("mdsquish.hist: Saved '%s'" % real_filename)                
    
    def network(self, filename=None, low_threshold=0.5, hi_threshold=0.9, cols=None, label_fontsize=8,
        max_links=9999999, labels=True, node_size=100, edges=True, save_gml=False, layout="neato",
        mark_clusters=False, cluster_alpha_back=0.8, cluster_node_size=3000, node_alpha=0.6, nodes=True,
        cluster_alpha_back2=1.0, mark_path=False, expected_branches=None, title=None, title_font_size=12,
        log_correct_node_sizes=None, bracket=None, edge_pad=0.03, 
        **kargs):
        """
        **Purpose**
            Draw a mdsquish network
        
        **Arguments**
            filename (Required)
                the filename to save the image to.
                               
            low_threshold (Optional, default=0.5)
                the correlation score will be used to build strength of links between the nodes (samples).
                
                However you will want to trim low strength connections between nodes otherwise all samples
                will end up connected to each other.
            
            hi_threshold (Optional, default=0.9)
                Thrshold to draw thick links
                
            max_links (Optional, default=9999999)
                put an arbitrary number of maximum links for each node as a limit. The algorithm always takes the five
                best nodes.
            
            cols (Optional, default=None)
                A list of colours to use for the nodes. 
            
            labels (Optional, default=True)
                draw labels?
                
            edges (Optional, default=True)
                draw the edges between nodes
                
            nodes (Optional, default=True)
                draw the nodes
            
            label_fontsize (Optional, default=7)
                The font size for the labels on the network
                
            node_size (Optional, default=100)
                The size of the spot for each node.
                
                mdsquish specific commang:
                You can also send a string for the node name and the mdsquish will collect the
                expression values and use them to plot the size of the nodes.
                
            bracket (Optional, default=False)
                Use a bracket value for the node_sizes. Will clamp the node size range to between this pair of values
                can help with the presentation, but may obscure larger changes. 
 
            log_correct_node_sizes (Optional, default=False)
                Often you may use log2 to transform the expression data, but this can give highly misleading node sizes
                Setting this to a log base will cause mdsquish to convert the value back into its
                unlogged value.
 
            node_cmap (Optional, default=False)
                mdsquish specific command:
                You can send a matplotlib colormap to use to color the nodes to mark intensity

            node_alpha (Optional, default=0.6)
                Alpha blending for the nodes.

            save_gml (Optional, default=False)
                A filename to save a GML file to. Suitable for loading into something like Cytoscape
                with the 'graphgmlreader' plugin

            layout (Optional, default="neato")
                A Graphviz layout engine. By default "neato" is used.
                
                Valid methods:
                
                "neato", "dot", "fdp", "sfdp", "circo"
                
                Note that not all Graphviz layouts are available, see also network.rooted()

            mark_path (Optional, default=False)
                mark paths on the network based on some method:
                
                    'longest_path'
                        Find the longest possible path between two nodes on the network.
                        (Actually the longest shortest path between any two nodes).
                        
                    'most_travelled' NOT IMPLEMENTED
                        The longest most travelled path in the network (Basically the sum of
                        all paths)
                    
                    'branches' 
                        Use all of the branches. Also expects a 'expected_branches' keyword argument
                    
                    'minimum_spanning_tree' NOT IMPLEMENTED
                        See minimum_spanning_tree() from networkx, this is the subgraph of the network
                        for a weighted tree. Uses Kruskal's algorithm. Sort of like longest_path
                        but it takes into account the weights.
                
            expected_branches (Optional, default=None)
                If mark_path is set to 'branches' then you must also set this variable to 
                the number of expected branches.
                
            mark_clusters (Optional, default=False)
                Divide the samples into an expected number of clusters.
                To enable this feature send an expected number of clusters in the sample.
                See also: mark_clusters()   
            
            cluster_alpha_back (Optional, default=0.8)
                Alpha color for the background groups if mark_clusters is used.
                
            cluster_alpha_back2 (Optional, default=1.0)
                Alpha color of the cluster background nodes.
                Semi undocumented feature, use cluster_alpha_back in preference.
            
            cluster_node_size (Optionale, default=3000)
                The node sizes of the cluster groups.
                
            title (Optional, default=None)
                An optional title to add to the plot.
            
            edge_pad (Optional, default=0.03)
                Fraction to pad around the edge of the network
                        
            draw_node_boundary (Optional, default=False)
                draw a boundary of nodes on the network. The nodes must be specifed by either node_coundary
                or you must have some sort of mark_path 
                
            node_boundary (Optional, default=None)
                a list of nodes to mark on the boundary.
        """
        assert self.data, "mdsquish.network(): run squish() first"
        if cols: assert len(cols) == len(self.names), "'cols' is not the same length as the number of conditions"
        if mark_path: assert mark_path in self.path_func_mapper, "network(): mark_path mode '%s' not found" % mark_path
        
        if (not self.__layout_data) or self.__last_thresholds != (low_threshold, hi_threshold, max_links):# Must regenerate if any of these change
            self.G = nx.Graph()
            for cind, row in enumerate(self.distances):
                self.G.add_node(self.names[cind], color=cols[cind])
        
            for cind, row in enumerate(self.distances):           
                scores = list(zip(self.names, row)) # The highest will always be self.
            
                scores.sort(key=operator.itemgetter(1), reverse=True)
                
                for item in scores:
                    rind, t = (self.names.index(item[0]), item[1])
                    if rind != cind: # stop self edges:
                        if t >= low_threshold:
                            if not self.G.has_edge(self.names[rind], self.names[cind]):
                                this_links = len(self.G.neighbors(self.names[rind]))
                                that_links = len(self.G.neighbors(self.names[cind]))
                                if (this_links < max_links and that_links < max_links): # or (this_links < min_links or that_links < min_links): 
                                    self.G.add_edge(self.names[rind], self.names[cind], weight=self.traversal_weight+(1.0-t))   # This is really for pathfinding. Penalty self.traversal_weight for moving each 1 cell.
        
            self.__last_thresholds = (low_threshold, hi_threshold, max_links)
            self.__layout_data =  nx.drawing.nx_agraph.graphviz_layout(self.G, layout) # Bug/Feature in NX 1.11
        else:
            config.log.info('mdsquish.conditions: Used preexisting network')
        
        # This is getting reused, parck into network_support helper?
        # Potentially could add above, but is a little tricky to maintain and probably fairly fast:
        # I need to repack node_size as mdsquish supports just passing the gene name as a string:
        node_color = "grey"
        if isinstance(node_size, str):
            ee = self.parent.find(node_size)
            assert ee, "'%s' not found" % node_size
            assert ee["conditions"], "no expression data found"

            if title is None or title is False: # If you want to remove the title by sending an empty string the empty string will pass an if not title: test. 
                title = node_size
            node_size = ee["conditions"]
        
            if log_correct_node_sizes:
                # correct the condition data for log
                node_size = [log_correct_node_sizes**(v) for v in node_size]
        
            node_size_t = numpy.array(node_size)
            if bracket:
                node_size_t[node_size_t>bracket[1]] = bracket[1]
                node_size_t -= bracket[0]
                node_size_t[node_size_t<0] = 0
                ma = bracket[1] - bracket[0]
                mi = 0
            else:
                ma = max(node_size)
                mi = min(node_size)
            
            # If no bracket then it will basically do a row_norm. With the bracket, this will scale it into 
            # A size range of 0 ... 2000
            # Scale the node sizes to make them aesthetically pleasing
            #print node_size_t
            node_size_t = ((node_size_t - mi) / (ma - mi)) * 2000 # puts it into range 0..1
            #print node_size_t
            node_size_t = [int(i) for i in node_size_t]
        
            # I need to reorder the node_sizes to match the nodes.
            # I don't know why this is required here, where commonly it is not required...
            ss = dict(list(zip(self.parent.getConditionNames(), node_size_t)))
            node_size = []
            for node in self.G.nodes():
                node_size.append(ss[node])
        
            if "node_cmap" in kargs and kargs["node_cmap"]:
                cmap = cm.get_cmap(kargs["node_cmap"], 2001) 
                cmap = cmap(numpy.arange(2001)) 
                cols = [cmap[i] for i in node_size]   
        
        return_data = network_support.unified_network_drawer(self.G, self.distances, self.names, filename=filename, 
            low_threshold=low_threshold, hi_threshold=hi_threshold, cols=cols, label_fontsize=label_fontsize,
            max_links=max_links, labels=labels, node_size=node_size, edges=edges, save_gml=save_gml, layout=layout,
            mark_clusters=mark_clusters, cluster_alpha_back=cluster_alpha_back, cluster_node_size=cluster_node_size,
            node_alpha=node_alpha, nodes=nodes, cluster_alpha_back2=cluster_alpha_back2, mark_path=mark_path,
            expected_branches=expected_branches, title=title, title_font_size=title_font_size, 
            traversal_weight=self.traversal_weight, width_adjuster=5,
            edge_pad=edge_pad, layout_data=self.__layout_data,
            **kargs)
        
        config.log.info("mdsquish.network: saved '%s'" % return_data["actual_filename"])
        
        return return_data
        
    def find_correlated_genes(self, mode="longest_path", R=0.7, patterns=None, **kargs):
        """
        **Purpose**
            find genes that correlate with a selection of patterns 
            across some path in the network as specified by 'mode'
            
        **Arguments**
            mode (Optional, default=False)
                Use a particular algorithm to find a certain type of path across the network.
                Note that 'weak' links are also considered for building the path
                
                Valid modes are:
                    'longest_path'
                        Find the longest possible path between two nodes on the network.
                        (Actually the longest shortest path between any two nodes).
                        
                    'most_travelled' NOT IMPLEMENTED
                        The longest most travelled path in the network (Basically the sum of
                        all paths)
                    
                    'branches' 
                        Use all of the branches. Also expects a 'expected_branches' keyword argument
                    
                    'minimum_spanning_tree' NOT IMPLEMENTED
                        See minimum_spanning_tree() from networkx, this is the subgraph of the network
                        for a weighted tree. Uses Kruskal's algorithm. Sort of like longest_path
                        but it takes into account the weights.

            R (Optional, default=0.7)
                The Pearson R score to consider a correlation as a match
                               
        **Returns**
            A dictionary in the form:
            {"path_1": glbase.expression, 
                "path_2": glbase.expression, 
                ... 
                "path_n": glbase.expression}

        """            
        assert mode in self.path_func_mapper, "find_correlated_genes(): mode '%s' not found" % mode
        assert self.G, "find_correlated_genes(): You will need to run mdsquish.network() first to populate the network"
        if mode == "branches":
            assert 'expected_branches' in kargs, "find_correlated_genes(): 'expected_branches' argument required when mode is branches"
        
        if not patterns:
            patterns = ["up", "dn", "updnup", "dnupdn"]
            
        nodes, edges = self.path_func_mapper[mode](self.G, **kargs) # get the path I'm going to use
                    
        results = {}
        for path_number, path_nodes in enumerate(nodes):
            #print path_number, path_nodes
            this_path_results = {"up": [], "dn": []}
            path_idxs = [self.parent.getConditionNames().index(i) for i in path_nodes]
        
            for gene in self.parent:
                # could also add updnup and dnupdn
                centroids = {"up": numpy.arange(len(path_nodes)), # needs to be regenerated depending upon the size of the path
                    "dn": numpy.arange(len(path_nodes))[::-1]}
                    
                for cen in centroids:
                    expns_on_path = [gene["conditions"][idx] for idx in path_idxs]
                    for expn in [expns_on_path, numpy.log2(numpy.array(expns_on_path)+0.1), numpy.log10(numpy.array(expns_on_path)+0.1)]:
                    
                        #print expn, gene["name"]
                        if True in [i > 0 for i in expns_on_path]: # only do if >0 in at least one entry
                            score = scipy.stats.pearsonr(centroids[cen], expn)[0]
                        
                            if score > R: 
                                this_path_results[cen].append(gene)   
                                #print score, centroids[cen], expn
                                break # only add once
            results["path_%s" % path_number] = this_path_results
            
        return results

    def get_path(self, start_node, end_node, mode='dijkstra', neighbours=0, filename=None, 
        low_threshold=0.5, hi_threshold=0.9, cols=None, label_fontsize=8,
        max_links=9999999, labels=True, node_size=100, edges=True, save_gml=False, layout="neato",
        mark_clusters=False, cluster_alpha_back=0.8, cluster_node_size=3000, node_alpha=0.6, nodes=True,
        cluster_alpha_back2=1.0, mark_path=False, expected_branches=None, title=None, title_font_size=12,
        log_correct_node_sizes=None, bracket=None, edge_pad=0.03, 
        **kargs):
        """
        **Purpose**
            Collect all of the samples from start_node to end_node by
            the most direct route (shortest path). Will take weight into consideration
        
        **Arguments**
            start_node (Required)
                the name of the startingnode
            
            end_node (Required)
                the name of the endingnode

            neighbours (Optional, default=0)
                collect n highest scoring neighbours to each node on the path. Will remove
                duplicates.
                
            filename (Optional, default=None)
                If true save a network marking the path.
                Will also support all arguments from mdsquish.network()
                
            mode (Optional, default='dijkstra')
                The path finding algorithm to use.
        """
        assert start_node, "mdsquish.get_path: Must specifiy start_node"
        assert end_node, "mdsquish.get_path: Must specifiy end_node"
        assert self.G, 'mdsquish.get_path: Need to run network() first to build a network'
               
        if mode == 'astar':
            def heuristic(u, v):
                return 0
        
            direct_path_names = nx.astar_path(self.G, start_node, end_node, heuristic=heuristic, weight='weight')  
        elif mode == 'dijkstra':
            #direct_path_names = nx.dijkstra_path(self.G, start_node, end_node, weight='weight')  
            direct_path_names = nx.dijkstra_path(self.G, start_node, target=end_node, weight='weight')#[1][end_node]
        
        if neighbours:
            direct_path_names, direct_path = network_support.populate_path_neighbours(self.G, direct_path_names)
        else:
            direct_path = list(zip(direct_path_names, direct_path_names[1:])) #
        #print direct_path
        
        if filename:
            node_color = "grey" # stored in node attributes, will override this
                   
            if "node_cmap" in kargs and kargs["node_cmap"]:
                cmap = cm.get_cmap(kargs["node_cmap"], 2001) 
                cmap = cmap(numpy.arange(2001)) 
                cols = [cmap[i] for i in node_size]   
        
            return_data = network_support.unified_network_drawer(self.G, self.distances, self.names, filename=filename, 
                low_threshold=low_threshold, hi_threshold=hi_threshold, cols=cols, label_fontsize=label_fontsize,
                max_links=max_links, labels=labels, node_size=node_size, edges=edges, save_gml=save_gml, layout=layout,
                mark_clusters=mark_clusters, cluster_alpha_back=cluster_alpha_back, cluster_node_size=cluster_node_size,
                node_alpha=node_alpha, nodes=nodes, cluster_alpha_back2=cluster_alpha_back2, mark_path=direct_path,
                expected_branches=expected_branches, title=title, title_font_size=title_font_size, 
                edge_pad=edge_pad, layout_data=self.__layout_data, width_adjuster=5,
                traversal_weight=self.traversal_weight,
                **kargs)
                
            config.log.info("mdsquish.get_path: saved '%s'" % return_data["actual_filename"])
        
        return(direct_path_names)

    def get_sum_of_paths(self, start_node_label, end_node_label, exclude=None, mode='dijkstra', 
        neighbours=0, filename=None, 
        low_threshold=0.5, hi_threshold=0.9, cols=None, label_fontsize=8,
        max_links=9999999, labels=True, node_size=100, edges=True, save_gml=False, layout="neato",
        mark_clusters=False, cluster_alpha_back=0.8, cluster_node_size=3000, node_alpha=0.6, nodes=True,
        cluster_alpha_back2=1.0, mark_path=False, expected_branches=None, title=None, title_font_size=12,
        log_correct_node_sizes=None, bracket=None, edge_pad=0.03, must_appear_in_at_least_n_percent_of_paths=20,
        **kargs):
        """
        **Purpose**
            Collect all of the samples from start_node to end_node by
            the most direct route (shortest path). Will take weight into consideration.
            
        
        **Arguments**
            start_nodes (Required)
                a label to use to collect starting nodes
            
            end_nodes (Required)
                a label to use to collect ending nodes
                
            exclude (Optional) 
                Specifically exclude these named nodes.

            must_appear_in_at_least_n_percent_of_paths (Optinoal, default=20)
                The node must appear in at least this percent paths to be counted.

            neighbours (Optional, default=0)
                collect n highest scoring neighbours to each node on the path. Will remove
                duplicates.
                
            filename (Optional, default=None)
                If true save a network marking the path.
                Will also support all arguments from mdsquish.network()
                
        """
        assert start_node_label, "mdsquish.get_sum_of_paths: Must specifiy start_node_label"
        assert end_node_label, "mdsquish.get_sum_of_paths: Must specifiy end_node_label"
        assert self.G, 'mdsquish.get_sum_of_paths: Need to run network() first to build a network'

        exclude = [] if not exclude else set(exclude)
        start_nodes = []
        end_nodes = []
        for c in self.G.nodes():
            if c in exclude:
                continue
            if start_node_label in c:
                start_nodes.append(c)
            elif end_node_label in c:
                end_nodes.append(c)

        assert start_nodes, "mdsquish.get_sum_of_paths: no samples found with start_node_label (%s)" % start_node_label
        assert end_nodes, "mdsquish.get_sum_of_paths: no samples found with end_nodes (%s)" % end_node_label

        # This is actually pretty fast...
        least_weighted_paths = []
        p = progressbar(len(start_nodes))
        for i, s in enumerate(start_nodes):
            for e in end_nodes:
                # get the least weighted path between these two nodes
                least_weighted_paths.append(nx.dijkstra_path(self.G, s, target=e, weight='weight'))
            p.update(i)

        # Build a new network from the subset of the old.
        # get all nodes that appear in at least t% of 
        must_appear_in_at_least_n_percent_of_paths = must_appear_in_at_least_n_percent_of_paths/100.0
        unwound = [item for sublist in least_weighted_paths for item in sublist] # flatten a 2D list to 1D           
        counted = collections.Counter(unwound)
        set_of_nodes_for_tree = []
        max_paths = len(end_nodes) * len(start_nodes)
        min_paths = min([len(end_nodes), len(start_nodes)])
        for k in counted:
            perc = (counted[k] - min_paths) / float(max_paths)
            #print k, perc, counted[k]
            if perc >= must_appear_in_at_least_n_percent_of_paths: # i.e. in more than 1 path
                set_of_nodes_for_tree.append(k)
        #print counted
        set_of_nodes_for_tree = set(set_of_nodes_for_tree)

        if neighbours:
            set_of_nodes_for_tree, direct_path = network_support.populate_path_neighbours(self.G, set_of_nodes_for_tree, degree=neighbours)

        # I have the nodes, but I need to put them back into order.
        # Cut out the appropriate part of the main network and then
        # get the longest, most direct path across the subnetwork
        newG = nx.Graph()
        # put in all the nodes
        for node in set_of_nodes_for_tree:
            newG.add_node(node)
        # Get out all the edges
        for node1 in set_of_nodes_for_tree:
            for node2 in set_of_nodes_for_tree:
                if (
                    node1 != node2
                    and self.G.has_edge(node1, node2)
                    and not newG.has_edge(node1, node2)
                ):
                    newG.add_edge(node1, node2, self.G.get_edge_data(node1, node2))
        #print newG.nodes()
        #print newG.edges()
        #mst = nx.minimum_spanning_tree(newG, weight='weight')
        # Get the final best path:
        best_path_nodes, best_path_edges = network_support.longest_path(newG) # Djikstra
        if neighbours: # Add neighbours back in if being used
            best_path_nodes, best_path_edges = network_support.populate_path_neighbours(self.G, best_path_nodes, degree=neighbours)

        #tree_path = mst.edges() # [i[0] for i in mst.edges()] + [mst.edges()[-1][1]]

        #print 'Path:', best_path_nodes, best_path_edges

        if filename:
            node_color = "grey" # stored in node attributes, will override this

            if "node_cmap" in kargs and kargs["node_cmap"]:
                cmap = cm.get_cmap(kargs["node_cmap"], 2001) 
                cmap = cmap(numpy.arange(2001)) 
                cols = [cmap[i] for i in node_size]   

            return_data = network_support.unified_network_drawer(self.G, self.distances, self.names, filename=filename, 
                low_threshold=low_threshold, hi_threshold=hi_threshold, cols=cols, label_fontsize=label_fontsize,
                max_links=max_links, labels=labels, node_size=node_size, edges=edges, save_gml=save_gml, layout=layout,
                mark_clusters=mark_clusters, cluster_alpha_back=cluster_alpha_back, cluster_node_size=cluster_node_size,
                node_alpha=node_alpha, nodes=nodes, cluster_alpha_back2=cluster_alpha_back2, mark_path=best_path_edges,
                traversal_weight=self.traversal_weight, width_adjuster=5,
                expected_branches=expected_branches, title=title, title_font_size=title_font_size, 
                edge_pad=edge_pad, layout_data=self.__layout_data, 
                **kargs)

            config.log.info("mdsquish.get_sum_of_paths: saved '%s'" % return_data["actual_filename"])

        return(best_path_nodes)

    def heatmap(self, filename=None, **kargs):
        """
        **Purpose**
            draw a heatmap of the Euclidean distances
            
        **Arguments**
            filename (Required)
                filename to save heatmap to.
                
            Will also accept almost all common heatmap arguments (see expression.heatmap())
        """
        assert self.data, "mdsquish.heatmap(): run squish() first"
        assert filename, "mdsquish.heatmap(): no filename!"

        if 'bracket' not in kargs:
            kargs['bracket'] = [0, 1]

        results = self.parent.draw.heatmap(filename=filename, data=self.distances, square=True,
            aspect='square', row_names=self.names, col_names=self.names, 
            colbar_label="Euclidean distance", **kargs)
            
        config.log.info("mdsquish.heatmap(): Save '%s'" % filename)
        return(None)
        
    def scatter(self, x, y, filename=None, spot_cols=None, label=False, alpha=0.8, 
        spot_size=40, label_font_size=7, cut=None, squish_scales=False, **kargs):
        """
        **Purpose**
            plot a scatter plot of cond1 against cond2.
        
        **Arguments**
            x, y (Required)
                PC dimension to plot as scatter
                Note that PC begin at 1 (and not at zero, as might be expected)
            
            filename (Required)
        
            spot_cols (Optional, default="black" or self.set_cols())
                list of colours for the samples, should be the same length as 
                the number of conditions. 
            
            label (Optional, default=False)
                label each spot with the name of the condition
                
            alpha (Optional, default=0.8)
                alpha value to use to blend the individual points
                
            spot_size (Optional, default=40)
                Size of the spots on the scatter
                
            label_font_size (Optional, default=7)
                Size of the spot label text, only valid if label=True
        
            cut (Optional, default=None)
                Send a rectangle of the form [topleftx, toplefty, bottomrightx, bottomrighty], cut out all of the items within that
                area and return their label and PC score    
                
            squish_scales (Optional, default=False)
                set the limits very aggressively to [minmin(x), minmax(y)]
        
        **Returns**
            None
            You can get PC data from pca.get_uvd()
        """
        assert filename, "scatter: Must provide a filename"     

        ret_data = None
        xdata = self.__transform[x-1]
        ydata = self.__transform[y-1]

        if "aspect" not in kargs:
            kargs["aspect"] = "square"

        fig = self.parent.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        cols = 'grey'
        if spot_cols:
            cols = spot_cols            

        ax.scatter(xdata, ydata, s=spot_size, alpha=alpha, edgecolors="none", c=cols)
        if label:
            for i, lab in enumerate(self.names):
                ax.text(xdata[i], ydata[i], lab, size=label_font_size, ha="center", va="top")

        # Tighten the axis
        if squish_scales:
            if "xlims" not in kargs:
                ax.set_xlim([min(xdata), max(xdata)])

            if "ylims" not in kargs:
                ax.set_ylim([min(ydata), max(ydata)])

        ax.set_xlabel("PC%s" % (x,)) # can be overridden via do_common_args()
        ax.set_ylabel("PC%s" % (y,))

        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])

        if cut:
            rect = matplotlib.patches.Rectangle(cut[0:2], cut[2]-cut[0], cut[3]-cut[1], ec="none", alpha=0.2, fc="orange")
            ax.add_patch(rect)

            tdata = []
            for i in range(len(xdata)):
                if (
                    xdata[i] > cut[0]
                    and xdata[i] < cut[2]
                    and ydata[i] < cut[1]
                    and ydata[i] > cut[3]
                ):
                    if self.rowwise: # grab the full entry from the parent genelist
                        dat = {"pcx": xdata[i], "pcy": ydata[i]}
                        dat.update(self.parent.linearData[i])
                        tdata.append(dat)
                    else:
                        tdata.append({"name": label[i], "pcx": xdata[i], "pcy": ydata[i]})
            if tdata:
                ret_data = genelist()
                ret_data.load_list(tdata)

        self.parent.draw.do_common_args(ax, **kargs)

        real_filename = self.parent.draw.savefigure(fig, filename)
        config.log.info("scatter: Saved 'PC%s' vs 'PC%s' scatter to '%s'" % (x, y, real_filename))
        return(ret_data)   