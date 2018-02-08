# parameter setting for almost all plotting
import matplotlib 

class Setparm:

    def __init__(self, parms):
        self.parmdict = {}
        for key in parms.keys():
            for i in range(0,len(parms[key])-1):
                self.parmdict[parms[key][i]] = parms[key][-1]	

    def setrc(self):
        matplotlib.rcParams.update(self.parmdict)
