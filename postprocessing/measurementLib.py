import numpy as np
import pdb
import matplotlib.pyplot as plt




class stretchedMeasurement(object):

    def __init__(self, cellEvents, subEvents, cellFreq, subFreq, cellNaturalFreq, dt, Delta, title):

        self.dt = dt
        self.cellEvents = cellEvents
        self.subEvents = subEvents
        self.cellFreq = cellFreq
        self.subFreq = subFreq
        self.cellNaturalFreq = cellNaturalFreq
        self.Delta = Delta
        self.title = title
        
        # Preprocessing:
        self.genTime()

        self.cellIx = self.genIndices(self.cellEvents)
        self.subIx = self.genIndices(self.subEvents)

        self.cellTheta = self.phaseGen(self.cellIx, self.t)
        self.subTheta = self.phaseGen(self.subIx, self.t)

        #self.clipTime()
        self.relativePhase()
        self.relativePhaseVelocityVsRelativePhase()

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

    def relativePhaseVelocity(self):
        self.dTheta2_dt = np.diff(self.dTheta2)/self.dt

        for i in range(0,self.dTheta2_dt.size):
            if self.dTheta2_dt[i]*self.dt > 0.5:
                self.dTheta2_dt[i] = (self.dTheta2[i+1] - self.dTheta2[i] - 1)/self.dt
            elif self.dTheta2_dt[i]*self.dt < - 0.5:
                self.dTheta2_dt[i] = (self.dTheta2[i+1] - self.dTheta2[i] - 1)/self.dt

    #######################################################
        
    def relativePhaseVelocityVsRelativePhase(self):

        self.relativePhaseVelocity()
        
        nBins  = 100

        dTheta2 = self.dTheta2[:-1]

        bins = np.linspace(0,1,nBins+1)
        dThetaRange = np.linspace(0,1,nBins+1)
        
        
        
        dTheta2_dt_accum = []
        for i in range(0,nBins):
            dTheta2_dt_accum.append([])

        for i in range(0,nBins):
            set1 = np.zeros_like(dTheta2)
            set2 = np.zeros_like(dTheta2)
            set1[np.where(dTheta2>bins[i])]=1
            set2[np.where(dTheta2<=bins[i+1])]=1
            dTheta2_dt_accum[i]= self.dTheta2_dt[np.where(set1*set2==1)[0]]
                

        self.dTheta_dt_of_dTheta = np.zeros(len(dTheta2_dt_accum))
        for i in range(0, len(dTheta2_dt_accum)):
            if len(dTheta2_dt_accum[i]) == 0:
               self.dTheta_dt_of_dTheta[i] = 0
            else:
               self.dTheta_dt_of_dTheta[i] = np.mean(dTheta2_dt_accum[i])

        self.dTheta_dt_phase = bins[0:-1]

        plt.plot(self.dTheta_dt_phase, self.dTheta_dt_of_dTheta)
        plt.show()
            




        


    #####################################################

    def hist_substratePhasesOnContraction(self):
        # Use this to get a histogram of substrate phases on initiation of cell contraction

        n, bins, patches = plt.hist(self.subTheta[self.cellIx])
        return n, bins, patches


    #####################################################

    def plotRelativePhase(self, figsize=(4.5,4.5), top=0.9, bottom=0.14, left=0.14, right=0.93):

        fig = plt.figure(figsize = figsize)
        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right)
        ax1 = fig.add_subplot(111)
        
        ax1.plot(self.t2, self.dTheta2)#, t2, dtheta2, t3, dtheta3, t4, dtheta4)
    
        ax1.set_ylim([0,1])
        ax1.set_yticks(np.linspace(0,1,3))
        ax1.set_ylabel(r"$\phi_{cell} - \phi_{sub}$", fontsize=16)
        ax1.set_xlabel("t (sec)")
        ax1.set_title(r"$\Delta$=" + '%.2f' % self.Delta, fontsize=16)

        plt.savefig(self.title + 'relativePhase.eps', dpi=160, facecolor='w')
        plt.show()
    
    
    

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


        if subTheta.size > self.cellTheta.size:
            subTheta = subTheta[0:self.cellTheta.size]
        elif subTheta.size < self.cellTheta.size:
            self.cellTheta = self.cellTheta[0:subTheta.size]

        dTheta = np.mod(self.cellTheta - subTheta, 1)

        minIndex = int(np.max([np.min(self.cellIx), np.min(subIx)]))
        maxIndex = int(np.min([np.max(self.cellIx), np.max(subIx)]))

        dTheta2 = dTheta[minIndex:maxIndex]
        t2 = self.t[minIndex:maxIndex]

        return t2, dTheta2













