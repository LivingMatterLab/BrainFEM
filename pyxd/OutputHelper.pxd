from ModelContainer cimport *
from DataContainer cimport *
from Element cimport *

cpdef ParaviewOutput(ModelContainer mc, DataContainer dc, int it = *)

cpdef PlotOutput(ModelContainer mc, DataContainer dc)

cpdef PlotFace(DataContainer dc, Element el, list n, int ndim)

cpdef WriteToLog(ModelContainer mc, str text)

cpdef WriteToOutput(ModelContainer mc, str filename, str text)
