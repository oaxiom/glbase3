"""
history.py

container for historical tracking of the lists.

"""

import sys, os, time

class historyItem:
    def __init__(self, parent, message):
        """
        """
        self.message = message
        self.time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        self.oldhash = parent.__hash__()

    def __str__(self):
        return("%s:%20s" % (self.time, self.message))

class historyContainer:
    def __init__(self, parent):
        self.l = []
        self.parent = parent

    def append(self, message):
        new_item = historyItem(self.parent, message)
        self.l.append(new_item)

    def __str__(self):
        if self.l:
            out = ["History for %s" % self.parent.name] + ["--------------------------"] + ["::%s" % item for item in self.l]
            out = "\r\n".join(out)
            return(out)
        else:
            return("History is Empty")
