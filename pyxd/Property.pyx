# -*- coding: utf-8 -*-
cdef class Property(object):
    def __init__(self,propertyType):
        self.localID = -1;
        self.type = propertyType;
        
    def __str__(self):
        return "Property " + str(self.localID) +":  \t type = " + self.type
