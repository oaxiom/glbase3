"""

network_support

Shared code for mdsquish and network, part of glbase

"""

import operator

import numpy, scipy
from matplotlib.collections import LineCollection
from matplotlib.colors import colorConverter, Colormap
import matplotlib.pyplot as plot
import matplotlib.cm as cm
import matplotlib.cbook as cb
import matplotlib.patches
import networkx as nx
from networkx.utils import is_string_like
from operator import itemgetter

from . import config, utils
from .draw import draw

gldraw = draw() # make the glbase draw system avaiablae

def draw_nodes(G, pos, ax=None, nodelist=None, node_size=300, node_col_override=None, node_color=None, node_shape='o',
    alpha=1.0, cmap=None, vmin=None, vmax=None, linewidths=None, label=None, zorder=2, **kargs):
    """
    Taken from draw_networkx_nodes() and modified for zorder, and a lot of boilerplate removed or simplified.
    
    Heavily modified these days though...
    """
    if nodelist is None: # set it to use all nodes if nodelist is None. is, as sometimes you get a numpy array
        nodelist = [n for n in G.nodes(data=True)] # Strip the networkx view business
        #nodelist = list(G) # data?
    elif isinstance(nodelist, list): # The node_boundary just sends back a list of node names
        # Convert to a tuple-like list of nodes, to match the output from G.nodes()
        nodelist = [(n, G.node[n]) for n in nodelist] # get the node back out from the full network    
    
    # set the colors from the attributes if present:
    if 'color' in nodelist[0][1]: # Test a node to see if color attrib present
        node_color = []
        for n in nodelist:
            node_color.append(n[1]['color'])
            #print n, n[1]['color']
    
    # override the colors set in the node attribs.
    if node_col_override:
        node_color = node_col_override
    
    # Fill in some colors if got to here and still no colors
    if not node_color:
        #print pos
        node_color = ["grey"] * len(nodelist)
    
    # Populate size from 'size attribute' if present.
    if 'size' in nodelist[0][1]:
        node_size = []
        for n in nodelist:
            node_size.append((n[1]['size'],))
    #else: Assume they know what they are doing
    
    xy = numpy.asarray([pos[v[0]] for v in nodelist])
    
    node_collection = ax.scatter(xy[:,0], xy[:,1], 
        s=node_size, 
        c=node_color,
        marker=node_shape, 
        cmap=cmap, 
        vmin=vmin, vmax=vmax, 
        alpha=alpha,
        linewidths=linewidths, 
        edgecolors=None,
        label=label)

    #plot.sci(node_collection)
    node_collection.set_zorder(zorder) # I need this modification to change the ordering

    return node_collection

def draw_edges(G, pos, ax, edgelist=None, width=1.0, width_adjuster=50, edge_color='k', style='solid',
    alpha=None, edge_cmap=None, edge_vmin=None, edge_vmax=None, traversal_weight=1.0, 
    edge_delengthify=0.15, 
    arrows=True,label=None, zorder=1, **kwds):
    """
    Code cleaned-up version of networkx.draw_networkx_edges
    
    New args:
    
    width_adjuster - the line width is generated from the weight if present, use this adjuster to thicken the lines (multiply)
    """
    if edgelist is None:
        edgelist = G.edges()

    if not edgelist or len(edgelist) == 0:  # no edges!
        return None

    # set edge positions
    edge_pos = [(pos[e[0]], pos[e[1]]) for e in edgelist]
    new_ep = []
    for e in edge_pos:
        x, y = e[0]
        dx, dy = e[1]
        
        # Get edge length
        elx = (dx - x) * edge_delengthify
        ely = (dy - y) * edge_delengthify
        
        x += elx
        y += ely
        dx -= elx
        dy -= ely
    
        new_ep.append(((x, y), (dx,dy)))
    edge_pos = numpy.asarray(new_ep)

    if not cb.iterable(width):
        #print [G.get_edge_data(n[0], n[1])['weight'] for n in edgelist]
        # see if I can find an edge attribute:
        if 'weight' in G.get_edge_data(edgelist[0][0], edgelist[0][1]): # Test an edge
            lw = [0.5+((G.get_edge_data(n[0], n[1])['weight']-traversal_weight)*width_adjuster) for n in edgelist]
        else:
            lw = (width,)
    else:   
        lw = width

    if not is_string_like(edge_color) and cb.iterable(edge_color) and len(edge_color) == len(edge_pos):
        if numpy.alltrue([cb.is_string_like(c) for c in edge_color]):
            # (should check ALL elements)
            # list of color letters such as ['k','r','k',...]
            edge_colors = tuple([colorConverter.to_rgba(c, alpha) for c in edge_color])
        elif numpy.alltrue([not cb.is_string_like(c) for c in edge_color]):
            # If color specs are given as (rgb) or (rgba) tuples, we're OK
            if numpy.alltrue([cb.iterable(c) and len(c) in (3, 4) for c in edge_color]):
                edge_colors = tuple(edge_color)
            else:
                # numbers (which are going to be mapped with a colormap)
                edge_colors = None
        else:
            raise ValueError('edge_color must consist of either color names or numbers')
    else:
        if is_string_like(edge_color) or len(edge_color) == 1:
            edge_colors = (colorConverter.to_rgba(edge_color, alpha), )
        else:
            raise ValueError('edge_color must be a single color or list of exactly m colors where m is the number or edges')

    edge_collection = LineCollection(edge_pos, 
        colors=edge_colors,
        linewidths=lw, 
        antialiaseds=(1,), 
        linestyle=style,
        transOffset=ax.transData, 
        zorder=zorder)
 
    edge_collection.set_label(label)
    ax.add_collection(edge_collection)

    if cb.is_numlike(alpha):
        edge_collection.set_alpha(alpha)

    if edge_colors is None:
        if edge_cmap is not None:
            assert(isinstance(edge_cmap, Colormap))
        edge_collection.set_array(numpy.asarray(edge_color))
        edge_collection.set_cmap(edge_cmap)
        if edge_vmin is not None or edge_vmax is not None:
            edge_collection.set_clim(edge_vmin, edge_vmax)
        else:
            edge_collection.autoscale()

    # update view
    '''
    minx = numpy.amin(numpy.ravel(edge_pos[:,:,0]))
    maxx = numpy.amax(numpy.ravel(edge_pos[:,:,0]))
    miny = numpy.amin(numpy.ravel(edge_pos[:,:,1]))
    maxy = numpy.amax(numpy.ravel(edge_pos[:,:,1]))

    w = maxx-minx
    h = maxy-miny
    padx,  pady = 0.05*w, 0.05*h
    corners = (minx-padx, miny-pady), (maxx+padx, maxy+pady)
    ax.update_datalim(corners)
    ax.autoscale_view()
    '''
    return(edge_collection)

def draw_node_labels(G, pos, labels=None, font_size=12, font_color='k',
    font_family='sans-serif', font_weight='normal', alpha=1.0, bbox=None, ax=None,
    zorder=1, **kargs):
    """
    **Purpose**
        Bug fix in networks.draw_networkx_labels - does not respect zorder
        
        ax is now required
    """
    assert ax, 'draw_node_labels: You must specify an axis to plot on'
    
    if labels is None:
        labels = dict((n, n) for n in G.nodes())

    # set optional alignment
    horizontalalignment = kargs.get('horizontalalignment', 'center')
    verticalalignment = kargs.get('verticalalignment', 'center')

    text_items = {}  # there is no text collection so we'll fake one
    for n, label in list(labels.items()):
        x, y = pos[n]
        #if not cb.is_string_like(label): # Assume users are a bit more savvy
        label = str(label)  # this will cause "1" and 1 to be labeled the same
            
        t = ax.text(x, y, label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            transform=ax.transData,
            bbox=bbox,
            clip_on=True,
            zorder=zorder) # bug fix
        text_items[n] = t

    return text_items

def hierarchical_clusters(G, data_table, node_names, expected_group_number):
    """
    **Purpose**
        Return a list of node names divided into the number of expected clusters.
        
        Clusters are generated by hierarchical clustering, then choosing the branch with the
        desired number of groups
        
    **Arguments**
        G (Required)
            The network
    
        data_table (Required)
            A table of distances
    
        node_names (Required)
            the list of node names (in order)
    
        expected_group_number (Required)
            You will need to estimate the number of groups in the data yourself.
            The algorithm will then dived the samples into the number of desired groups
        
    **Returns**
        A dictionary of lists, in the form:
        
        {"cluster_1": [sam1, sam2 ...], "cluster_2": [sam1, sam2 ...] ... "cluster_n"}
        
    """
    assert G, "hierarchical_clusters: No network specified"
    
    xs = numpy.arange(0.1, 1.1, 0.001)[::-1]
    
    dist = scipy.spatial.distance.pdist(data_table) 
    linkage = scipy.cluster.hierarchy.linkage(dist, method="complete")

    for threshold in xs: # move up the threshold till I find the desired number of groups
        # dist.max() is WRONG! Y also has a minimum...
        clus = scipy.cluster.hierarchy.fcluster(linkage, threshold*dist.max(), 'distance') # Make sure this matches in do_network()
        num_componets = len(set(clus))
        if num_componets >= expected_group_number:
            break
            
    # It is possible to reach here with the wrong number of groups, although probably extremely rare.
    if num_componets != expected_group_number:
        config.log.warning("mark_clusters(): Could not find '%s' number of groups in data, using closest last group number '%s'" % (expected_group_number, num_componets))
    
    # map the clusters back to names:
    clusters = {}
    for idx, cluster in enumerate(clus):
        clus_name = "cluster_%s" % cluster
        if clus_name not in clusters:
            clusters[clus_name] = []
        clusters[clus_name].append(node_names[idx])
        
    return(clusters)
    
def longest_path(G, **kargs):
    """
    **Purpose**
        return a list of node names and edgelists for the longest, most direct path
        in the network
        
    **Arguments** 
        None
        
    **Returns**
    """
    assert G, "longest_path: No network specified"
    
    longest_path = []
    for node1 in G:
        for node2 in G:
            try:
                n = nx.dijkstra_path(G, node1, node2, weight='weight') # i.e. nx.shortest_path()
            except nx.exception.NetworkXNoPath: # Fail silently
                n = [] # i.e. don't add to longest path
            if len(n) > len(longest_path):
                longest_path = n

    # convert nodes to edgelist:
    longest_edges = list(zip(longest_path, longest_path[1:]))   
    return(longest_path, longest_edges) # does not support multiple paths

def __path_similarity(p1, p2):
    """
    test a path for % similarity
    
    Very rough and ready, just measures the number of nodes in common, not the order or edges.
    
    I should check the Diez paper for a better arrangement.
    """
    A = set(p1)
    B = set(p2)
    AB = A & B
    
    perc_sim = len(AB) / float(min(len(A),len(B))) * 100.0
    return(perc_sim)    

def branches(G, **kargs):
    """
    **Purpose**
        return a list of the the expected number of branches in this network
        
    **Arguments** 
        expected_branches
            The estimated number of expected branches to recover
        
    **Returns**
    """
    assert G, "branches: No network specified"
    
    longest_nodes, longest_edges = longest_path(G) # first get the longest path
    longest_nodes = longest_nodes[0] # unpack
    longest_edges = longest_edges[0]
        
    longest_branches = []
    for origin_node in longest_nodes:
        # find all of the longest paths from this node:
        for node2 in G:  
            try:
                n = None
                n = nx.shortest_path(G, origin_node, node2)        
            except nx.exception.NetworkXNoPath:
                pass
            # work out if the branch crosses back onto the longest_path, and trim the path at the point it rejoins the
            # longest_path (this stops paths doubling along the natural optimal path) 
            if n:
                sofar = []
                for r in n[1:]: # Dont include the first node!
                    if r in longest_nodes:
                        break # bug out of path
                    else:
                        sofar.append(r)
                if len(sofar) > 2: # trim trivial branches
                    longest_branches.append(sofar)
                
    # remove duplicate paths and append a length key:
    longest_branches = list(set([tuple(l) for l in longest_branches]))
    longest_branches = [(len(l), l) for l in longest_branches]
        
    # sort out the longest branches:
    longest_branches.sort(key=operator.itemgetter(0), reverse=True)
    longest_branches = [l[1] for l in longest_branches] # kill the score:

    # remove lower-scoring subset paths
    survivors = []
    worst_score = 0.0
    for i1, p1 in enumerate(longest_branches):
        for i2, p2 in enumerate(longest_branches):
            if i2 < i1: # only check items above this one
                worst_score = max(worst_score, __path_similarity(p1, p2))
        if worst_score < 51.0:
            survivors.append(p1)
        worst_score = 0.0    
    longest_branches = survivors
    
    # Now select the expected_branches from the longest_branches
    ret_nodes = [longest_nodes] + longest_branches[0:kargs["expected_branches"]]
    ret_edges = [list(zip(n, n[1:])) for n in ret_nodes]

    #longest_edges = zip(longest_path, longest_path[1:])   
    return(ret_nodes, ret_edges)

def minimum_spanning_tree(G, **kargs):
    """
    **Purpose**
        Find the minimum weighted spanning tree across the network
        
        See networkx:     minimum_spanning_tree
    
    **Arguments**
        None
        
    """
    assert G, "minimum_spanning_tree(): No network specified"
    nodes = nx.minimum_spanning_tree(self.G) # returns a new network. Will need to unpack
    nodes = nodes.nodes()
    edges = list(zip(nodes, nodes[1:])) # This is wrong?
    return([nodes], [edges])

path_func_mapper = {"longest_path": longest_path,
    "minimum_spanning_tree": minimum_spanning_tree,
    "branches": branches
    }
    
def unified_network_drawer(G, correlation_table, names, filename=None, low_threshold=0.5, hi_threshold=0.9, 
    cols=None, label_fontsize=8, edge_alpha=1.0, trim_isolated_nodes=False,
    max_links=9999999, labels=True, node_size=100, edges=True, save_gml=False, layout="neato",
    mark_clusters=False, cluster_alpha_back=0.8, cluster_node_size=3000, node_alpha=0.6, nodes=True,
    cluster_alpha_back2=1.0, mark_path=None, mark_paths=None, path_color='red', title=None, edge_pad=0.03, title_font_size=12,
    traversal_weight=0.0, draw_node_boundary=False, node_boundary=None,
    width_adjuster=20, # default for MDSquish 
    layout_data=None, # preexisting layout data
    **kargs): 
    """
    In the kargs:
    expected_branches
    edge_color
    edge_width
    
    unified network draw system.
    
    In use by:
    network.conditions()
    network.genes()
    mdsquish.network()
    
    TODO: 
    The edge_width is a bit messy. 
    There are three major arguments:
    edge_width, width_adjuster, traversal_weight 
    and they interact in complicated ways.
    
    zorder, lowest is further back higher is further forward
    
    """
    # Kargs and defaults:
    edge_color = 'grey'
    edge_width = 1.0
    if 'edge_color' in kargs and kargs['edge_color']: 
        edge_color = kargs['edge_color']
    if 'edge_width' in kargs and kargs['edge_width']: 
        edge_width = kargs['edge_width']
    
    # optional return data
    ret_groups = None
    ret_nodes = None
    ret_edges = None

    if layout_data: 
        pos = layout_data
    else:
        pos = nx.drawing.nx_agraph.graphviz_layout(G, layout) # Bug in NX 1.11
    
    # trim isolated nodes
    if trim_isolated_nodes:
        # The problem is, all the attributes are also unsynced...
        pass
        
    fig = gldraw.getfigure(**kargs)
    ax = fig.add_subplot(111)

    # get cols back in the node order:
    # Nice Py2.7 line put back to uglier style.
    #sam_map = {cond: cols[idx] for idx, cond in enumerate(self.getConditionNames())} # reorder
    if cols:
        sam_map = dict((cond, cols[idx]) for idx, cond in enumerate(names))
        cols = [sam_map[cond] for cond in G]
    else:
        cols = "grey"
    
    if node_boundary: # Make the background nodes not in the node boundary more transparent
        node_alpha = 0.1
    
    if nodes:
        draw_nodes(G, pos, ax=ax, 
        node_size=node_size, 
        node_color=cols, 
        alpha=node_alpha, 
        linewidths=0, zorder=5
        )
    
    #print 'univerted:', [2.0-(i[2]['weight']) for i in G.edges(data=True)]
    elarge = [(u,v,d) for (u,v,d) in G.edges(data=True) if ((traversal_weight+1.0)-d['weight']) >= hi_threshold] # I pad and invert the weight so that pathfinding works correctly
    esmall = [(u,v,d) for (u,v,d) in G.edges(data=True) if ((traversal_weight+1.0)-d['weight']) < hi_threshold] # valid as all edges must be less than 1.0-hi

    # mark clusters
    if mark_clusters:
        groups = hierarchical_clusters(G, correlation_table, names, mark_clusters)  
        
        # get a colormap for the groups:
        colormap = cm.get_cmap("Set3", len(groups)+1) 
        colormap = colormap(numpy.arange(len(groups)+1)) 

        # draw the groups by size?
        gsizes = {g: len(groups[g]) for g in groups} # Assume py2.7 now...
        
        for g in groups:
            node_color = utils.rgba_to_hex(colormap[int(g.replace("cluster_", ""))-1])
            draw_nodes(G, pos, ax=ax, nodelist=groups[g], 
                node_size=cluster_node_size, 
                node_col_override=node_color, 
                node_color=node_color, 
                alpha=cluster_alpha_back2, 
                linewidths=0, 
                zorder=-gsizes[g])
        # Draw an alpha box over the entire network to fade out the groups
        # This could be replaced by imshow for nicer effect.
        xl = ax.get_xlim()
        yl = ax.get_ylim()
        if cluster_alpha_back:
            ax.add_patch(matplotlib.patches.Rectangle((xl[0], yl[0]), xl[1]-xl[0], yl[1]-yl[0], facecolor="white", edgecolor='none', zorder=0, alpha=cluster_alpha_back))
        ret_groups = groups
        
    # edges
    if edges: 
        draw_edges(G, pos, ax, edgelist=elarge, width=edge_width, width_adjuster=width_adjuster, alpha=edge_alpha, edge_color='#666666', traversal_weight=traversal_weight, zodrder=4)
        draw_edges(G, pos, ax, edgelist=esmall, width=edge_width, width_adjuster=width_adjuster, alpha=edge_alpha/2.0, edge_color='#bbbbbb', traversal_weight=traversal_weight, zorder=3)

    # labels
    if labels:
        draw_node_labels(G, pos, ax=ax, font_size=label_fontsize, font_family='sans-serif', zorder=5)

    if mark_path:       
        if isinstance(mark_path, list): # ou are probably sending your own path
            draw_edges(G, pos, ax, edgelist=mark_path, width=5.0, alpha=1.0, edge_color=path_color, 
                width_adjuster=width_adjuster*2.0, traversal_weight=traversal_weight, zorder=6) # in front of nodes
        else:
            ret_nodes, ret_edges = path_func_mapper[mark_path](G, **kargs) # call the appropriate function
            cmap = cm.get_cmap("gist_ncar", len(ret_edges)) 
            cmap = cmap(numpy.arange(len(ret_edges))) 
            color = [utils.rgba_to_hex(cmap[i]) for i, e in enumerate(ret_edges)]
            
            for i, e in enumerate(ret_edges):
                draw_edges(G, pos, ax, edgelist=e, width=3.0, alpha=1.0, edge_color=color[i], traversal_weight=traversal_weight, zorder=6)
                # Don't draw the nodes, you also would need to get out the node name and properly reorder the node_size
            mark_path = ret_nodes # For compatibility with network_boundary
    elif mark_paths:
        for p in mark_paths:
            draw_edges(G, pos, ax, edgelist=mark_path, width=5.0, alpha=0.5, edge_color=path_color, 
                width_adjuster=300, traversal_weight=traversal_weight, zorder=6) # in front of nodes
    
    # node_boundary
    if draw_node_boundary:
        if node_boundary:
            boundary = node_boundary # I assume it is already a boundary
        elif mark_path: # use a path if no network_boundary sent
            boundary = nx.node_boundary(G, mark_path) 
        else:
            raise AssertionError('asked to draw a boundary, but no network_boundary or path available')
        draw_nodes(G, pos, ax=ax, 
            nodelist=boundary,
            node_size=node_size*1.2, #don't use node_color, see if the draw_nodes can pick it up from the attributes
            alpha=0.9, linewidths=0.0, zorder=3)  

    # clean up matplotlib gubbins:
    ax.set_position([0,0,1,1])
    ax.set_frame_on(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    # make nice edges (by default it chooses far too generous borders):
    xy = numpy.asarray([pos[v] for v in G.nodes()])
    x_min, x_max = min(xy[:,0]), max(xy[:,0])
    y_min, y_max = min(xy[:,1]), max(xy[:,1])
    x_pad = (x_max - x_min) * edge_pad
    y_pad = (y_max - y_min) * edge_pad
    ax.set_xlim(x_min-x_pad, x_max+x_pad)
    ax.set_ylim(y_min-y_pad, y_max+y_pad)
    
    if title:
        #ax.set_title(title)
        ax.text(x_min-(x_pad//2), y_min-(y_pad//2), title, ha='left', size=title_font_size)
           
    if save_gml:
        nx.write_gml(G, save_gml)
        config.log.info("network_drawer: saved GML '%s'" % save_gml)
    
    actual_filename = gldraw.savefigure(fig, filename)
    
    # Load the return data:
    ret = {"actual_filename": actual_filename}
    if ret_groups:
        ret["groups"] = ret_groups
    if mark_path:
        ret["nodes"] = ret_nodes
        ret["edges"] = ret_edges
    
    return(ret)

def populate_path_neighbours(G, path_names, degree=1):
    """
    get all of the n degree nearest neighbours and repopulate the path list.
    
    returns a tuple containing ([nodelist], [edgelist])
    
    """
    last_node = None
    new_path = []
    node_pool = []
    for node_name in path_names:
        # I need to add all of the nodes on the path first, to stop excessive pruning of neighbours
        node_pool.append(node_name)
    node_pool = set(node_pool) # So I can look up if it already on the path and so need to look for a better neighbour
    
    new_direct_path_names = []
    for node_name in path_names:    
        if last_node:
            # add the main path
            new_path.append((last_node, node_name))
            new_direct_path_names.append(node_name)
        # get all neighbours
        this_node = G[node_name]
        # sorted by weight:
        ss = []
        for k in this_node:
            ss.append((k, this_node[k]['weight']))
        ss = sorted(ss, key=itemgetter(1)) # Lower is now better
        
        neighbours_done = 0
        for n in ss:
            if n[0] not in node_pool:
                new_path.append((node_name, n[0]))
                node_pool.add(n[0])
                new_direct_path_names.append(n[0])
                neighbours_done += 1
                if neighbours_done >= degree:
                    break
        
        last_node = node_name
    direct_path_names = new_direct_path_names # repack for return data
    direct_path = new_path # Already in edge format.
    return(direct_path_names, direct_path)