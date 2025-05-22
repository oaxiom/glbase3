"""

network

Only loaded if networkx is available.

TODO:
. genes() and conditions() started, see network_support.unified_network_drawer
. BUG: max_links does no use the best links first. I think this is fixed?
"""

import operator

import numpy, heapq, itertools
import matplotlib
import matplotlib.pyplot as plot
import matplotlib.cm as cm
import networkx as nx # Should only be possible to get here if networkx is available.

from . import config
from . import network_support
from . import utils

class network:
    def __init__(self, parent):
        """
        network class, designed not to be called directly, but like this:

        expn.network.genes()

        although I suppose you could do this:

        expn = expression()
        netw = network(expn)

        netw.genes(...)

        If you really wanted to.

        **Argumants**
            parent must be an expression object.

        """
        self.parent = parent

    def __get_network_size(self, mode, names=None):
        if mode == "conditions":
            length = len(self.parent.getConditionNames())
        elif mode == "genes":
            length = len(self.parent[names])

        return length

    def __corr_network(self, lo, hi, length, max_links, mode="conditions", names="name", cols=None):
        """
        Helper function for building the correlation network

        returns the network object and the correlation_table
        """
        if cols and isinstance(cols, list):
            cols = cols
        elif cols:
            cols = [cols] * length
        else:
            cols = ["grey"] * length

        # build the network;
        G = nx.Graph()

        if mode == "conditions":
            correlation_table = numpy.corrcoef(self.parent.getExpressionTable().T)[:,:-1]
            #correlation_table *= correlation_table # i.e. R^2

            self.names = self.parent.getConditionNames()
            for cind, row in enumerate(correlation_table):
                if cols:
                    G.add_node(self.names[cind], color=cols[cind]) # each sample
                else:
                    G.add_node(self.names[cind])

        elif mode == "genes":
            correlation_table = numpy.corrcoef(self.parent.getExpressionTable())
            #correlation_table *= correlation_table # i.e. R^2

            # Names cannot be empty, or they come back as none. So I hack in "-" to stand for unlabbelled
            self.names = self.parent[names]
            for i, n in enumerate(self.names):
                if n == "":
                    self.names[i] = "-"

            for cind, row in enumerate(correlation_table):
                G.add_node(self.names[cind], color=cols[cind]) # each gene

        # Build the links:
        for cind, row in enumerate(correlation_table):
            scores = list(zip(self.names, row)) # The highest will always be self.

            scores.sort(key=operator.itemgetter(1), reverse=True)
            #print scores

            done = 0
            for item in scores:
                rind, t = (self.names.index(item[0]), item[1])
                if rind != cind: # stop self edges:
                    if t >= lo:
                        if not G.has_edge(self.names[rind], self.names[cind]):
                            if len(list(G.neighbors(self.names[rind]))) < max_links and len(list(G.neighbors(self.names[cind]))) < max_links:
                                G.add_edge(self.names[rind], self.names[cind], weight=1.0-t)
                                done += 1

                            if done > max_links:
                                break

        self.correlation_table = correlation_table

        return self.names, G, correlation_table, cols

    def __trim_isolated_nodes(self, G, cols, node_size):
        """
        Remove isolated nodes from the network
        """
        deg = G.degree()
        if not isinstance(cols, int):
            node_size = list(node_size) # Just in case you numpy'd it

        to_remove = [n for n in deg if deg[n] == 0]
        to_remove_idx = [n[0] for n in enumerate(deg) if deg[n[1]] == 0]
        G.remove_nodes_from(to_remove)
        to_remove_idx.sort()
        to_remove_idx.reverse()
        for idx in to_remove_idx: # safe if done from the highest to smallest
            if isinstance(cols, list):
                del cols[idx]
            if isinstance(cols, list):
                del node_size[idx]
        return G, cols, node_size

    def genes(self, filename=None, cols=None, label_fontsize=8,
        edge_alpha=1.0,
        max_links=9999999, labels=True, node_size=100, edges=True, name_key="name", save_gml=False, layout="neato",
        trim_isolated_nodes=False,
        mark_clusters=False, cluster_alpha_back=0.8, cluster_node_size=3000, node_alpha=0.6, nodes=True,
        cluster_alpha_back2=1.0, log_correct_node_sizes=False,
        low_threshold=None, hi_threshold=None,
        **kargs):
        """
        **Purpose**
            Plot a network based on the (R^2) correlation coefficient based on co-expreesion similariy of genes.

            Although you can do this method on the entire gene expression, generally figures with more than
            a few hundred nodes will be impossible to interpret

            NOTE: This method requires networkx and pygraphviz are both available.

        **Arguments**
            filename (Required)
                the filename to save the image to.

            low_threshold (Required)
                the correlation score will be used to build strength of links between the nodes (samples).

                However you will want to trim low strength connections between nodes otherwise all samples
                will end up connected to each other.

            hi_threshold (Required)
                Thrshold to draw thick links

            max_links (Optional, default=9999999)
                put an arbitrary number of maximum links for each node as a limit. The algorithm takes the max_links
                best nodes.

            cols (Optional, default=None)
                A list of colours to use for the nodes.

            node_size (Optional, default=None)
                A list of node sizes. or just a single node size for all nodes

            labels (Optional, default=True)
                draw labels?

            edges (Optional, default=True)

            edge_width (Optional, default=1.5)
                the width of the edges, in points.

            edge_alpha (Optional, default=1.0)
                the transparency of the edges, lo_threshold edges will be drawn at
                edge_alpha/2.0

            label_fontsize (Optional, default=8)
                The font size for the labels on the network

            name_key (Optional, default="name")
                Label to use for the nodes

            trim_isolated_nodes (Optional, default=False)
                remove nodes with no edges

            save_gml (Optional, default=False)
                A filename to save a GML file to. Suitable for loading into something like Cytoscape
                with the 'graphgmlreader' plugin

            layout (Optional, default="neato")
                A Graphviz layout engine. By default "neato" is used.

                Valid methods:

                "neato", "dot", "fdp", "sfdp", "circo"

                Note that not all Graphviz layouts are available.

            mark_clusters (Optional, default=False)
                Divide the samples into an expected number of clusters.
                To enable this feature send an expected number of clusters in the sample.
                See also: mark_clusters()

        **Returns**
            None
        """
        assert low_threshold, "You must specifiy a low_threshold"
        assert hi_threshold,  "You must specifiy a hi_threshold"
        assert layout in ("neato", "dot", "fdp", "sfdp", "circo"), "%s layout not supported" % layout
        assert len(set(self.parent[name_key])) == len(self.parent[name_key]), "The node names are not unique"
        if isinstance(cols, list):
            assert len(self.parent[name_key]) == len(cols), "The number of nodes and the length of 'cols' do not match"

        length = self.__get_network_size(mode="genes", names=name_key)
        if length <= 1:
            config.log.warning("Cannot generate a network with %s items" % length)
            return None

        conds, G, correlation_table, cols = self.__corr_network(low_threshold, hi_threshold, length, max_links, mode="genes", names=name_key, cols=cols)

        return_data = network_support.unified_network_drawer(G, correlation_table, self.names, filename=filename,
            low_threshold=low_threshold, hi_threshold=hi_threshold, cols=cols, label_fontsize=label_fontsize,
            max_links=max_links, labels=labels, node_size=node_size, edges=edges, save_gml=save_gml, layout=layout,
            mark_clusters=mark_clusters, cluster_alpha_back=cluster_alpha_back, cluster_node_size=cluster_node_size,
            node_alpha=node_alpha, nodes=nodes, cluster_alpha_back2=cluster_alpha_back2, width_adjuster=2.0,
            **kargs)

        config.log.info("network.conditions(): saved '%s'" % return_data["actual_filename"])

        return return_data

    def conditions(self, filename=None, cols=None, label_fontsize=8,
        max_links=9999999, labels=True, node_size=100, edges=True, save_gml=False, layout="neato",
        mark_clusters=False, cluster_alpha_back=0.8, cluster_node_size=3000, node_alpha=0.6, nodes=True,
        cluster_alpha_back2=1.0, log_correct_node_sizes=False, bracket=None,
        low_threshold=None, hi_threshold=None, title=None, title_font_size=12,
        **kargs):
        """
        **Purpose**
            Plot a network based on the (R^2) correlation coefficient for all pairs of
            samples in the expression object.

            NOTE: This method requires networkx and pygraphviz are both available.

        **Arguments**
            filename (Required)
                the filename to save the image to.

            low_threshold (Required)
                the correlation score will be used to build strength of links between the nodes (samples).

                However you will want to trim low strength connections between nodes otherwise all samples
                will end up connected to each other.

            hi_threshold (Required)
                Thrshold to draw thick links

            max_links (Optional, default=9999999)
                put an arbitrary number of maximum links for each node as a limit. The algorithm always takes the five
                best nodes.

            cols (Optional, default=None)
                A list of colours to use for the nodes.

            labels (Optional, default=True)
                draw labels on the nodes

            edges (Optional, default=True)
                draw the edges between nodes

            edge_color (Optional, default=grey)
                The edge color to use.

            edge_width (Optional, default=1.0)
                edge width in points.

            nodes (Optional, default=True)
                draw the nodes

            label_fontsize (Optional, default=7)
                The font size for the labels on the network

            title (Optional, default=False)
                Print the title in the lower left corner.

            title_font_size (Optional, default=12)
                The font size for the title. In pts

            node_size (Optional, default=100)
                The size of the spot for each node.

                You can also send a string for the node name and the mdsquish will collect the
                expression values and use them to plot the size of the nodes.

            log_correct_node_sizes (Optional, default=False)
                Often you may use log2 to transform the expression data, but this can give highly misleading node sizes
                Setting this to a log base will cause netwok to convert the value back into its
                original unlogged value.

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

        **Returns**
            if mark_clusters is >0 then it returns the members for each of the groups.

            Otherwise None
        """
        assert low_threshold, "You must specifiy a low_threshold"
        assert hi_threshold,  "You must specifiy a hi_threshold"
        assert layout in ("neato", "dot", "fdp", "sfdp", "circo"), "%s layout not supported" % layout
        assert len(set(self.parent.getConditionNames())) == len(self.parent.getConditionNames()), "The node names are not unique"
        if isinstance(cols, list):
            assert len(self.parent.getConditionNames()) == len(cols), "The number of nodes and the length of 'cols' do not match"

        length = self.__get_network_size(mode="conditions")
        if length <= 1:
            config.log.warning("Cannot generate a network with %s items" % length)
            return None

        conds, G, correlation_table, cols = self.__corr_network(low_threshold, hi_threshold, length, max_links, cols=cols)
        # conds is redundant and should be removed; cols is redundant and should be removed

        # This is getting reused, parck into network_support helper?
        # I need to repack node_size as mdsquish supports just passing the gene name as a string:
        node_color = "grey"
        if isinstance(node_size, str):
            ee = self.parent.find(node_size)
            assert ee, "'%s' not found" % node_size
            assert ee["conditions"], "no expression data found"

            if not title:
                title = node_size
            node_size = ee["conditions"]

            if log_correct_node_sizes:
                # correct the condition data for log
                node_size = [log_correct_node_sizes**(v) for v in node_size]
                #print node_size

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
        elif isinstance(node_size, int):
            node_size_t = [node_size] * len(self.parent.getConditionNames())
        else:
            node_size_t = node_size

        # I need to reorder the node_sizes to match the nodes.
        # I don't know why this is required here, where commonly it is not required...
        ss = dict(list(zip(self.parent.getConditionNames(), node_size_t)))
        node_size = []
        for node in G.nodes():
            node_size.append(ss[node])

            if "node_cmap" in kargs and kargs["node_cmap"]:
                cmap = cm.get_cmap(kargs["node_cmap"], 2001)
                cmap = cmap(numpy.arange(2001))
                cols = [cmap[i] for i in node_size]

        return_data = network_support.unified_network_drawer(G, correlation_table, self.names, filename=filename,
            low_threshold=low_threshold, hi_threshold=hi_threshold, cols=cols, label_fontsize=label_fontsize,
            max_links=max_links, labels=labels, node_size=node_size, edges=edges, save_gml=save_gml, layout=layout,
            mark_clusters=mark_clusters, cluster_alpha_back=cluster_alpha_back, cluster_node_size=cluster_node_size,
            node_alpha=node_alpha, nodes=nodes, cluster_alpha_back2=cluster_alpha_back2, title=title,
            title_font_size=title_font_size, width_adjuster=2.0,
            **kargs)

        config.log.info("network.conditions: saved '%s'" % return_data["actual_filename"])

        return None

    def rooted(self, filename=None, threshold=0.5, max_links=2, cols=None, label_fontsize=8,
        root_sample=None, **kargs):
        """
        **Purpose**
            Plot a network based on the (R^2) correlation coefficient for all pairs of
            samples in the expression object.

            This variant plots a circular map with 'root_sample' at the centre and all other samples
            plotted in relation to that sample.

            NOTE: This method requires networkx is available.

            This should be folded back into network.conditions()?

        **Arguments**
            filename (Required)
                the filename to save the image to.

            threshold (Optional, default=0.5)
                The minimum threshold to consider a link.

            root_sample (Required)
                The sample to use as root.

            max_links (Optional, default=2)
                Take at most the top <num_links> to build the network.

            cols (Optional, default=None)
                A list of colours to use for the nodes.

            label_fontsize (Optional, default=7)
                The font size for the labels on the network

        **Returns**

        """
        assert root_sample, "correlation_rooted_network: You must specify a root_sample"
        assert root_sample in self.parent.getConditionNames(), "correlation_rooted_network: root_sample '%s' not found" % (root_sample)

        length = self.__get_network_size(mode="conditions")
        if length <= 1:
            config.log.warning("Cannot generate a network with %s items" % length)
            return None

        conds, G, correlation_table, cols = self.__corr_network(threshold, threshold, length, max_links, cols=cols)

        pos = nx.drawing.nx_agraph.graphviz_layout(G, prog='twopi', root="'%s'" % root_sample) # bug in names with spaces in it.

        fig = self.parent.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        # nodes
        # get cols back in the node order:
        #sam_map = {cond: cols[idx] for idx, cond in enumerate(self.getConditionNames())} # reorder
        if cols:
            sam_map = dict((cond, cols[idx]) for idx, cond in enumerate(self.parent.getConditionNames())) # put back in for py2.6 compat.
            cols = [sam_map[cond] for cond in G]
        else:
            cols = "red"

        nx.draw_networkx_nodes(G, pos, node_size=100, node_color=cols, alpha=0.6, linewidths=0)

        elarge = [(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] <= 1.0-threshold]

        # edges
        nx.draw_networkx_edges(G, pos, edgelist=elarge, width=1, alpha=0.5, edge_color='grey')

        # labels
        nx.draw_networkx_labels(G, pos, font_size=label_fontsize, font_family='sans-serif')

        # clean up matplotlib gubbins:
        ax.set_position([0,0,1,1])
        ax.set_frame_on(False)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        actual_filename = self.parent.draw.savefigure(fig, filename)

        config.log.info("correlation_rooted_network: saved '%s'" % actual_filename)
