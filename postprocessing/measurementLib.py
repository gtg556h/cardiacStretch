import numpy as np
import pdb
import matplotlib.pyplot as plt




class stretchedMeasurement(object):

    def __init__(self, cellEvents, subEvents, cellFreq, subFreq, cellNaturalFreq, dt, Delta):

        self.dt = dt
        self.cellEvents = cellEvents
        self.subEvents = subEvents
        self.cellFreq = cellFreq
        self.subFreq = subFreq
        self.cellNaturalFreq = cellNaturalFreq
        self.Delta = Delta
        
        # Preprocessing:
        self.genTime()

        self.cellIx = self.genIndices(self.cellEvents)
        self.subIx = self.genIndices(self.subEvents)

        self.cellTheta = self.phaseGen(self.cellIx, self.t)
        self.subTheta = self.phaseGen(self.subIx, self.t)

        #self.clipTime()
        self.relativePhase()


    ######################################################

    def genTime(self):
        self.t = np.arange(0, np.max([np.max(self.subEvents),np.max(self.cellEvents)]), self.dt)


    ######################################################
        
    def genIndices(self, events):

        ix = np.zeros(np.size(events))
        
        for i in range(0, np.size(events)):
            ix[i] = int(np.argmin(np.abs(events[i] - self.t)))

        return ix.astype(int)
        

    ######################################################
    
    def phaseGen(self,ix,t):
        phase = np.zeros(t.shape)
        ixDiff = np.diff(ix,n=1,axis=0)

        for ii in range(0,ix.shape[0]-1):
            #phase[ix[ii,0],0] = 0
            for jj in range(0,int(ixDiff[ii])+1):
                phase[int(ix[ii])+jj] = float(jj)/float(ixDiff[ii])

        return phase #, ixDiff


    #####################################################
    
    def clipTime(self):
        return 1


    #####################################################
    
    def relativePhase(self):
        self.dTheta = np.mod(self.cellTheta - self.subTheta, 1)

        self.minIndex = int(np.max([np.min(self.cellIx), np.min(self.subIx)]))
        self.maxIndex = int(np.min([np.max(self.cellIx), np.max(self.subIx)]))

        self.dTheta2 = self.dTheta[self.minIndex:self.maxIndex]
        self.t2 = self.t[self.minIndex:self.maxIndex]


    #####################################################

    def hist_substratePhasesOnContraction(self):
        # Use this to get a histogram of substrate phases on initiation of cell contraction

        n, bins, patches = plt.hist(self.subTheta[self.cellIx])
        return n, bins, patches


    #####################################################
    

class unstretchedMeasurement(object):

    def __init__(self, cellEvents, cellFreq, cellNaturalFreq, dt, startTime):

        self.dt = dt
        self.cellEvents = cellEvents
        self.cellFreq = cellFreq
        self.cellNaturalFreq = cellNaturalFreq
        self.startTime = startTime
        
        # Preprocessing:
        self.genTime()

        self.cellIx = self.genIndices(self.cellEvents)
        self.cellTheta = self.phaseGen(self.cellIx, self.t)

        
        #self.clipTime()
        #self.relativePhase()


    ######################################################

    def genTime(self):
        self.t = np.arange(0, np.max(self.cellEvents), self.dt)


    ######################################################

    def genIndices(self, events):

        ix = np.zeros(np.size(events))
        
        for i in range(0, np.size(events)):
            ix[i] = int(np.argmin(np.abs(events[i] - self.t)))

        return ix.astype(int)
        

    ######################################################
    
    def phaseGen(self,ix,t):
        phase = np.zeros(t.shape)
        ixDiff = np.diff(ix,n=1,axis=0)
        
        for ii in range(0,ix.shape[0]-1):
            #phase[ix[ii,0],0] = 0
            for jj in range(0,int(ixDiff[ii])+1):
                phase[int(ix[ii])+jj] = float(jj)/float(ixDiff[ii])

        return phase #, ixDiff


    def relativePhase(self, subTheta, subIx):
        dTheta = np.mod(self.cellTheta - subTheta, 1)

        minIndex = int(np.max([np.min(self.cellIx), np.min(subIx)]))
        maxIndex = int(np.min([np.max(self.cellIx), np.max(subIx)]))

        dTheta2 = dTheta[minIndex:maxIndex]
        t2 = self.t[minIndex:maxIndex]

        return dTheta2
