"""

markov.py

markov chain support for glbase

markov chains are made available throught this interface:

import markov.chain [as chain]

the markov chains here are specific to DNA but could probably be
expanded to further word-based items.

"""

from . import config

default_letter_set = ["a", "c", "g", "t"]

class _node:
    pass

class chain:
    def __init__(self, letter_set=default_letter_set):
        """
        **Purpose**
            Initialise a Markov chain object

        **Argumnets**
            letter_set (Required, default = [a,c,g,t])
                This is a hidden argument that defines all of the available
                letter sets
        """
        self.letter_set = letter_set
        self.__level = None

    def __str__(self):
        """
        (Internal)
        Return a string representation of the Markov Chain
        """

        return """
        Markov Chain:\n
        Levels: %s\n
        """ % (
            self.__len__(),
        )

    def __len__(self):
        """
        (Internal)
        Returns the markov chain level
        """
        if self.__level:
            return(self.__level)
        return(-1)
