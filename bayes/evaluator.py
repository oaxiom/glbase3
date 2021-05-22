"""Classes and functions for efficiently evaluating networks."""

from math import log
import random

import numpy as N

from . import data, cpd, prior, config, network
from .util import *

N.random.seed()

#
# Exceptions
#
class CyclicNetworkError(Exception):
    msg = "Network has cycle and is thus not a DAG."


#
# Localscore Cache
#
class LocalscoreCache(object):
    """ A LRU cache for local scores.

    Based on code from http://code.activestate.com/recipes/498245/
    """

    _params = (
        config.IntParameter(
            'localscore_cache.maxsize',
            "Max number of localscores to cache. Default=-1 means unlimited size.",
            default=-1
        )
    )

    def __init__(self, evaluator, cachesize=None):
        self._cache = {}
        self._queue = deque()
        self._refcount = {}
        self.cachesize = cachesize or config.get('localscore_cache.maxsize')

        self.neteval = evaluator
        self.hits = 0
        self.misses = 0

    def __call__(self, node, parents):
        _queue = self._queue
        _refcount = self._refcount
        _cache = self._cache
        _maxsize = self.cachesize

        index = tuple([node] +  parents)

        # get from cache or compute
        try:
            score = _cache[index]
            self.hits += 1
        except KeyError:
            score = _cache[index] = self.neteval._cpd(node, parents).loglikelihood()
            self.misses += 1

        # if using LRU cache (maxsize != -1)
        if _maxsize > 0:
            # record that key was accessed
            _queue.append(index)
            _refcount[index] = _refcount.get(index, 0) + 1

            # make variables local
            _len = len
            # purge LRU entry
            while _len(_cache) > _maxsize:
                k = _queue.popleft()
                _refcount[k] -= 1
                if not _refcount[k]:
                    del _cache[k]
                    del _refcount[k]

            # Periodically compact the queue by duplicate keys
            if _len(_queue) > _maxsize * 4:
                for _ in range(_len(_queue)):
                    k = _queue.popleft()
                    if _refcount[k] == 1:
                        _queue.append(k)
                    else:
                        _refcount[k] -= 1

        return score



#
# Network Evaluators
#
class NetworkEvaluator(object):
    """Base Class for all Network Evaluators.

    Provides methods for scoring networks but does not eliminate any redundant
    computation or cache retrievals.
    
    """

    def __init__(self, data_, network_, prior_=None, localscore_cache=None):

        self.network = network_
        self.data = data_
        self.prior = prior_ or prior.NullPrior()
        
        self.datavars = list(range(self.data.variables.size))
        self.score = None
        self._localscore = localscore_cache or LocalscoreCache(self)
        self.localscore_cache = self._localscore

    #
    # Private Interface
    # 
    def _globalscore(self, localscores):
        # log(P(M|D)) +  log(P(M)) == likelihood + prior
        return N.sum(localscores) + self.prior.loglikelihood(self.network)
    
    def _cpd(self, node, parents):
        #return cpd.MultinomialCPD(
            #self.data.subset(
                #[node] + parents,            
                #N.where(self.data.interventions[:,node] == False)[0])) 
        return cpd.MultinomialCPD(
            self.data._subset_ni_fast([node] + parents))


    def _score_network_core(self):
        # in this implementation, we score all nodes (even if that means
        # redundant computation)
        parents = self.network.edges.parents
        self.score = self._globalscore(
            self._localscore(n, parents(n)) for n in self.datavars
        )
        return self.score

    #
    # Public Interface
    #
    def score_network(self, net=None):
        """Score a network.

        If net is provided, scores that. Otherwise, score network previously
        set.

        """

        self.network = net or self.network
        return self._score_network_core()

    def alter_network(self, add=[], remove=[]):
        """Alter network by adding and removing sets of edges."""

        self.network.edges.add_many(add)
        self.network.edges.remove_many(remove)
        return self.score_network()
    
    def randomize_network(self): 
        """Randomize the network edges."""

        self.network = network.random_network(self.network.nodes)
        return self.score_network()

    def clear_network(self):     
        """Clear all edges from the network."""

        self.network.edges.clear()
        return self.score_network()


class SmartNetworkEvaluator(NetworkEvaluator):
    def __init__(self, data_, network_, prior_=None, localscore_cache=None):
        """Create a 'smart' network evaluator.

        This network evaluator eliminates redundant computation by keeping
        track of changes to network and only rescoring the changes. This
        requires that all changes to the network are done through this
        evaluator's methods. 

        The network can be altered by the following methods:
            * alter_network
            * score_network
            * randomize_network
            * clear_network

        The last change applied can be 'undone' with restore_network

        """

        super(SmartNetworkEvaluator, self).__init__(data_, network_, prior_, 
                                                    localscore_cache)

        # can't use this with missing data
        if self.data.missing.any():
            msg = "Cannot use the SmartNetworkEvaluator with missing data."
            raise Exception(msg)

        # these represent that state that we intelligently manage
        self.localscores = N.zeros((self.data.variables.size), dtype=float)
        self.dirtynodes = set(self.datavars)
        self.saved_state = None

    #
    # Private Interface
    #
    def _backup_state(self, added, removed):
        self.saved_state = (
            self.score,                     # saved score
            #[(n,self.localscores[n]) for n in self.dirtynodes],
            self.localscores.copy(),        # saved localscores
            added,                          # edges added
            removed                         # edges removed
        )

    def _restore_state(self):
        if self.saved_state:
            self.score, self.localscores, added, removed = self.saved_state
            #self.score, changedscores, added, removed = self.saved_state
        
        #for n,score in changedscores:
            #self.localscores[n] = score

        self.network.edges.add_many(removed)
        self.network.edges.remove_many(added)
        self.saved_state = None
        self.dirtynodes = set()

    def _score_network_core(self):
        # if no nodes are dirty, just return last score.
        if len(self.dirtynodes) == 0:
            return self.score

        # update localscore for dirtynodes, then re-calculate globalscore
        parents = self.network.edges.parents
        for node in self.dirtynodes:
            self.localscores[node] = self._localscore(node, parents(node))
        
        self.dirtynodes = set()
        self.score = self._globalscore(self.localscores)

        return self.score
    
    #
    # Public Interface
    #
    def score_network(self, net=None):
        """Score a network.

        If net is provided, scores that. Otherwise, score network previously
        set.

        """

        if net:
            add = [edge for edge in net.edges if edge not in self.network.edges]
            remove = [edge for edge in self.network.edges if edge not in net.edges]
        else:
            add = remove = []
        
        return self.alter_network(add, remove)
    
    def alter_network(self, add=[], remove=[]):
        """Alter the network while retaining the ability to *quickly* undo the changes."""

        # make the required changes
        # NOTE: remove existing edges *before* adding new ones. 
        #   if edge e is in `add`, `remove` and `self.network`, 
        #   it should exist in the new network. (the add and remove cancel out.
        self.network.edges.remove_many(remove)
        self.network.edges.add_many(add)    

        # check whether changes lead to valid DAG (raise error if they don't)
        affected_nodes = set(unzip(add, 1))
        if affected_nodes and not self.network.is_acyclic(affected_nodes):
            self.network.edges.remove_many(add)
            self.network.edges.add_many(remove)
            raise CyclicNetworkError()
        
        
        # accept changes: 
        #   1) determine dirtynodes
        #   2) backup state
        #   3) score network (but only rescore dirtynodes)
        self.dirtynodes.update(set(unzip(add+remove, 1)))
        self._backup_state(add, remove)
        self.score = self._score_network_core()

        return self.score
       
    def randomize_network(self):
        """Randomize the network edges."""

        newnet = network.random_network(self.network.nodes)
        return self.score_network(newnet)

    def clear_network(self):
        """Clear all edges from the network."""

        return self.alter_network(remove=list(self.network.edges))

    def restore_network(self):
        """Undo the last change to the network (and score).
        
        Undo the last change performed by any of these methods:
            * score_network
            * alter_network
            * randomize_network
            * clear_network

        """

        self._restore_state()
        return self.score


class GibbsSamplerState(object):
    """Represents the state of the Gibbs sampler.

    This state object can be used to resume the Gibbs sampler from a particaular point.
    Note that the state does not include the network or data and it's upto the caller to ensure
    that the Gibbs sampler is resumed with the same network and data.

    The following values are saved:
        - number of sampled scores (numscores)
        - average score (avgscore)
        - most recent value assignments for missing values (assignedvals)

    """

    def __init__(self, avgscore, numscores, assignedvals):
        self.avgscore = avgscore
        self.numscores = numscores
        self.assignedvals = assignedvals

    @property
    def scoresum(self):
        """Log sum of scores."""
        return self.avgscore + N.log(self.numscores)


class MissingDataNetworkEvaluator(NetworkEvaluator):
    #
    # Parameters
    #
    _params = (
        config.IntParameter(
            'gibbs.burnin',
            """Burn-in period for the gibbs sampler (specified as a multiple of
            the number of missing values)""",
            default=10
        ),
        config.StringParameter(
            'gibbs.max_iterations',
            """Stopping criteria for the gibbs sampler.
            
            The number of Gibb's sampler iterations to run. Should be a valid
            python expression using the variable n (number of missing values).
            Examples:

                * n**2  (for n-squared iterations)
                * 100   (for 100 iterations)
            """,
            default="n**2"
        )
    )

    def __init__(self, data_, network_, prior_=None, localscore_cache=None, 
                 **options): 
        """Create a network evaluator for use with missing values.

        This evaluator uses a Gibb's sampler for sampling over the space of
        possible completions for the missing values.

        For more information about Gibb's sampling, consult:

            1. http://en.wikipedia.org/wiki/Gibbs_sampling
            2. D. Heckerman. A Tutorial on Learning with Bayesian Networks. 
               Microsoft Technical Report MSR-TR-95-06, 1995. p.21-22.

       
        Any config param for 'gibbs' can be passed in via options.
        Use just the option part of the parameter name.

        """

        super(MissingDataNetworkEvaluator, self).__init__(data_, network_,
                                                         prior_)
        self._localscore = None  # no cache w/ missing data
        config.setparams(self, options)
        
    def _init_state(self):
        parents = self.network.edges.parents

        self.cpds = [self._cpd(n, parents(n)) for n in self.datavars]
        self.localscores = [cpd.loglikelihood() for cpd in self.cpds]
        self.data_dirtynodes = set(self.datavars)

    def _score_network_core(self):
        # update localscore for data_dirtynodes, then calculate globalscore.
        for n in self.data_dirtynodes:
            self.localscores[n] = self.cpds[n].loglikelihood()

        self.data_dirtynodes = set()
        self.score = self._globalscore(self.localscores)
        return self.score

    def _alter_data(self, row, col, value):
        oldrow = self.data.observations[row].copy()
        self.data.observations[row,col] = value

        # update data_dirtynodes
        affected_nodes = set(self.network.edges.children(col) + [col])
        self.data_dirtynodes.update(affected_nodes)

        # update cpds
        for node in affected_nodes:
            datacols = [node] + self.network.edges.parents(node)
            if not self.data.interventions[row,node]:
                self.cpds[node].replace_data(
                        oldrow[datacols],
                        self.data.observations[row][datacols])

    def _alter_data_and_score(self, row, col, value):
        self._alter_data(row, col, value)
        return self._score_network_core()

    def _calculate_score(self, chosenscores, gibbs_state):
        # discard the burnin period scores and average the rest
        burnin_period = self.burnin * \
                        self.data.missing[self.data.missing==True].size

        if gibbs_state:
            # resuming from a previous gibbs run. so, no burnin required.
            scoresum = logsum(N.concatenate((chosenscores, [gibbs_state.scoresum])))
            numscores = len(chosenscores) + gibbs_state.numscores
        elif len(chosenscores) > burnin_period:
            # remove scores from burnin period.
            nonburn_scores = chosenscores[burnin_period:]
            scoresum = logsum(nonburn_scores)
            numscores = len(nonburn_scores)
        else:
            # this occurs when gibbs iterations were less than burnin period.
            scoresum = chosenscores[-1]
            numscores = 1
        
        score = scoresum - log(numscores)
        return score, numscores

    def _assign_missingvals(self, indices, gibbs_state):
        if gibbs_state:
            assignedvals = gibbs_state.assignedvals
        else:
            arities = [v.arity for v in self.data.variables]
            assignedvals = [random.randint(0, arities[col]-1) for row,col in indices]
        
        self.data.observations[unzip(indices)] = assignedvals

    def score_network(self, net=None, gibbs_state=None):
        """Score a network.

        If net is provided, scores that. Otherwise, score network previously
        set.

        The default stopping criteria is to run for n**2 iterations.

        gibbs_state is the state of a previous run of the Gibb's sampler.  With
        this, one can do the following::
        
            myeval = evaluator.MissingDataNetworkEvaluator(...)
            myeval.score_network(...)
            gibbs_state = myeval.gibbs_state
            cPickle.dump(gibbs_state, 'gibbs_state.txt')

            # look at results, do other analysis, etc
            # If we decide that we need further Gibb's sampler iterations, we
            # don't need to restart
            gibbs_state = cPickle.load(open('gibbs_state.txt'))
            myeval = evaluator.MissingDataNetworkEvaluator(...)

            # continue with the previous run of the Gibb's sampler
            myeval.score_network(
                gibbs_state=gibbs_state,
                stopping_criteria=lambda i,N: i>200*N**2
            )

        """

        self.network = net or self.network

        # create some useful lists and local variables
        missing_indices = unzip(N.where(self.data.missing==True))
        num_missingvals = len(missing_indices)
        n = num_missingvals
        max_iterations = eval(self.max_iterations)
        arities = [v.arity for v in self.data.variables]
        chosenscores = []

        self._assign_missingvals(missing_indices, gibbs_state)
        self._init_state()

        # Gibbs Sampling: 
        # For each missing value:
        #    1) score net with each possible value (based on node's arity)
        #    2) using a probability wheel, sample a value from the possible values
        iters = 0
        while iters < max_iterations:
            for row,col in missing_indices:
                scores = [self._alter_data_and_score(row, col, val) \
                             for val in range(arities[col])]
                chosenval = logscale_probwheel(list(range(len(scores))), scores)
                self._alter_data(row, col, chosenval)
                chosenscores.append(scores[chosenval])
            
            iters += num_missingvals

        self.chosenscores = N.array(chosenscores)
        self.score, numscores = self._calculate_score(self.chosenscores, gibbs_state)

        # save state of gibbs sampler
        self.gibbs_state = GibbsSamplerState(
            avgscore=self.score, 
            numscores=numscores, 
            assignedvals=self.data.observations[unzip(missing_indices)].tolist()
        )

        return self.score


class MissingDataExactNetworkEvaluator(MissingDataNetworkEvaluator):
    """MissingDataNEtworkEvaluator that does an exact enumeration.

    This network evaluator enumerates over all possible completions of the
    missing values.  Since this is a combinatorial space, this class is only
    feasible with datasets with few missing values.

    """

    def score_network(self, net=None, stopping_criteria=None, gibbs_state=None):
        """Score a network.

        If net is provided, scores that. Otherwise, score network previously
        set.

        Note: See MissingDataNetworkEvaluator.score_network for more information
        about arguments.

        """

        self.network = net or self.network
        
        # create some useful lists and local variables
        missing_indices = unzip(N.where(self.data.missing==True))
        num_missingvals = len(missing_indices)
        possiblevals = [list(range(self.data.variables[col].arity)) for row,col in missing_indices]

        self._init_state()
        
        # Enumerate through all possible values for the missing data (using
        # the cartesian_product function) and score.
        scores = []
        for assignedvals in cartesian_product(possiblevals):
            for (row,col),val in zip(missing_indices, assignedvals):
                self._alter_data(row, col, val)
            scores.append(self._score_network_core())

        # average score (in log space)
        self.score = logsum(scores) - log(len(scores))
        return self.score


class MissingDataMaximumEntropyNetworkEvaluator(MissingDataNetworkEvaluator):
    """MissingDataNetworkEvaluator that uses a different space of completions.

    This evaluator only samples from missing value completions that result in a
    maximum entropy discretization for the variables with missing values. This
    is useful when the rest of the variables are maximum-entropy discretized
    because then all variables have the same entropy.

    """

    def _do_maximum_entropy_assignment(self, var):
        """Assign values to the missing values for this variable such that
        it has a maximum entropy discretization.
        """

        arity = self.data.variables[var].arity
        numsamples = self.data.samples.size

        missingvals = self.data.missing[:,var]
        missingsamples = N.where(missingvals == True)[0]
        observedsamples = N.where(missingvals == False)[0]
        
        # maximum entropy discretization for *all* samples for this variable
        numeach = numsamples/arity
        assignments = flatten([val]*numeach for val in range(arity))
        for i in range(numsamples - len(assignments)):  
            assignments.append(i)

        # remove the values of the observed samples
        for val in self.data.observations[observedsamples, var]:
            assignments.remove(val)

        N.random.shuffle(assignments)
        self.data.observations[missingsamples,var] = assignments


    def _assign_missingvals(self, missingvars, gibbs_state):
        if gibbs_state:
            assignedvals = gibbs_state.assignedvals
            self.data.observations[N.where(self.data.missing==True)] = assignedvals
        else:
            for var in missingvars:
                self._do_maximum_entropy_assignment(var)
 

    def _swap_data(self, var, sample1, choices_for_sample2):
        val1 = self.data.observations[sample1, var]

        # try swapping till we get a different value (but don't keep trying
        # forever)
        for _ in range(len(choices_for_sample2)/2):
            sample2 = random.choice(choices_for_sample2)
            val2 = self.data.observations[sample2, var]
            if val1 != val2:
                break

        self._alter_data(sample1, var, val2)
        self._alter_data(sample2, var, val1)

        return (sample1, var, val1, sample2, var, val2)
    
    def _undo_swap(self, row1, col1, val1, row2, col2, val2):
        self._alter_data(row1, col1, val1)
        self._alter_data(row2, col2, val2) 

    def score_network(self, net=None, gibbs_state=None):
        """Score a network.

        If net is provided, scores that. Otherwise, score network previously
        set.  
        
        Note: See MissingDataNetworkEvaluator.score_network for more information
        about arguments.

        """

        self.network = net or self.network

        # create some useful lists and counts
        num_missingvals = self.data.missing[self.data.missing == True].shape[0]
        n = num_missingvals
        max_iterations = eval(self.max_iterations)
        chosenscores = []
        
        # determine missing vars and samples
        missingvars = [v for v in self.datavars if self.data.missing[:,v].any()]
        missingsamples = [N.where(self.data.missing[:,v] == True)[0] \
                            for v in self.datavars]

        self._assign_missingvals(missingvars, gibbs_state)
        self._init_state()

        # iteratively swap data randomly amond samples of a var and score
        iters = 0
        while iters < max_iterations:
            for var in missingvars:  
                for sample in missingsamples[var]:
                    score0 = self._score_network_core()
                    swap = self._swap_data(var, sample, missingsamples[var]) 
                    score1 = self._score_network_core() 
                    chosenval = logscale_probwheel([0,1], [score0, score1])
                    
                    if chosenval == 0:
                        self._undo_swap(*swap)
                        chosenscores.append(score0)
                    else:
                        chosenscores.append(score1)

            iters += num_missingvals

        self.chosenscores = N.array(chosenscores)
        self.score, numscores = self._calculate_score(self.chosenscores, gibbs_state)
        
        # save state of gibbs sampler
        self.gibbs_state = GibbsSamplerState(
            avgscore=self.score, 
            numscores=numscores, 
            assignedvals=self.data.observations[
                N.where(self.data.missing==True)
            ].tolist()
        )
        
        return self.score

#
# Parameters
#
_pmissingdatahandler = config.StringParameter(
    'evaluator.missingdata_evaluator',
    """
    Evaluator to use for handling missing data. Choices include:
        * gibbs: Gibb's sampling
        * maxentropy_gibbs: Gibbs's sampling over all completions of the
          missing values that result in maximum entropy discretization for the
          variables.  
        * exact: exact enumeration of all possible missing values (only
                 useable when there are few missing values)
    """,
    config.oneof('gibbs', 'exact', 'maxentropy_gibbs'),
    default='gibbs'
)

_missingdata_evaluators = {
    'gibbs': MissingDataNetworkEvaluator,
    'exact': MissingDataExactNetworkEvaluator,
    'maxentropy_gibbs': MissingDataMaximumEntropyNetworkEvaluator
}

def fromconfig(data_=None, network_=None, prior_=None):
    """Create an evaluator based on configuration parameters.
    
    This function will return the correct evaluator based on the relevant
    configuration parameters.
    
    """

    data_ = data_ or data.fromconfig()
    network_ = network_ or network.fromdata(data_)
    prior_ = prior_ or prior.fromconfig()

    if not data_.missing.any():
        return SmartNetworkEvaluator(data_, network_, prior_)

    e = _missingdata_evaluators[config.get('evaluator.missingdata_evaluator')]
    return e(data_, network_, prior_)

