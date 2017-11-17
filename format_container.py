"""

The new format container, with docs!

"""

import copy

class fc():
    def __init__(self, name, format, description):
        self.data = format
        self.name = name
        
        if description:
            self.description = description
        
        elif "__description__" in self.data:
            self.description = self.data["__description__"]
        
        else:   
            self.description = "None"
        self.__doc__ = self.description
    
    def __str__(self):
        a = ["Name: %s" % self.name,
            "Format: %s" % self.data,
            "Description: %s" % self.description]
        return ("\n".join(a))
        
    def __getitem__(self, key):
        return(self.data[key])
    
    def __setitem__(self, key, value):
        self.data[key] = value
        
    def __repr__(self):
        return("<fc: %s" % self.name)
        
    def __iter__(self):
        for n in self.data:
            yield n
            
    def has_key(self, key):
        return(key in self.data)
        
    def __in__(self, key):
        return(key in self.data)
    
    def update(self, dictionary):
        self.data.update(dictionary)    

    def copy(self):
        c = copy.copy(self)
        c.data = self.data.copy() # dict copy
        return(c)