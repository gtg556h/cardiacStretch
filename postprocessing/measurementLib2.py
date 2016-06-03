import numpy as np
import pdb
import matplotlib.pyplot as plt



class measurement(object):

    def __init__(self, cellEvents, cellFreq, cellNaturalFreq, dt, startTime, title, maxTime):

        self.dt = dt
        self.cellEvents = cellEvents
        self.cellFreq = cellFreq
        self.cellNaturalFreq = cellNaturalFreq
        self.startTime = startTime
        self.title = title
        
        # Preprocessing:
        self.genTime(maxTime)

        self.cellIx = self.genIndices(self.cellEvents)
        self.cellTheta = self.phaseGen(self.cellIx, self.t)

    ######################################################

    def genTime(self, maxTime):
        #self.t = np.arange(0, np.max([np.max(self.subEvents),np.max(self.cellEvents)]), self.dt)
        self.t = np.arange(0, maxTime, self.dt)


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

    

##############################################################
##############################################################
##############################################################


class unstretchedMeasurement(measurement):

    def __init__(self, cellEvents, cellFreq, cellNaturalFreq, dt, startTime, title):

        maxTime = np.max(cellEvents)
        super().__init__(cellEvents, cellFreq, cellNaturalFreq, dt, startTime, title, maxTime)


    ######################################################
    
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



    
############################################################
############################################################
############################################################


class stretchedMeasurement(measurement):

    def __init__(self, cellEvents, cellFreq, cellNaturalFreq, dt, startTime, title, subEvents, subFreq, Delta):

        maxTime = np.max([np.max(cellEvents), np.max(subEvents)])
        super().__init__(cellEvents, cellFreq, cellNaturalFreq, dt, startTime, title, maxTime)

        self.subEvents = subEvents
        self.subFreq = subFreq
        self.Delta = Delta
        self.subIx = self.genIndices(self.subEvents)
        self.subTheta = self.phaseGen(self.subIx, self.t)
        self.relativePhase()

    #####################################################

    def relativePhase(self):
        self.dTheta = np.mod(self.cellTheta - self.subTheta, 1)

        self.minIndex = int(np.max([np.min(self.cellIx), np.min(self.subIx)]))
        self.maxIndex = int(np.min([np.max(self.cellIx), np.max(self.subIx)]))

        self.dTheta = self.dTheta[self.minIndex:self.maxIndex]
        self.t2 = self.t[self.minIndex:self.maxIndex]

        self.dTheta_t = self.relativePhaseVelocity(self.dTheta, self.dt)

        self.dTheta_t_bins, self.dTheta_t_dTheta = self.genRelativePhaseVelocityVsRelativePhase(self.dTheta, self.dTheta_t)

    #####################################################

    def relativePhaseVelocity(self, dTheta, dt):

        dTheta_t = np.diff(dTheta) / dt

        for i in range(0, dTheta_t.size):
            if dTheta_t[i] * dt > 0.5:
                dTheta_t[i] = (dTheta[i+1] - dTheta[i] - 1) / dt
            elif dTheta_t[i] * dt < -0.5:
                dTheta_t[i] = (dTheta[i+1] - dTheta[i] + 1) / dt

        return dTheta_t

    #######################################################

    def genRelativePhaseVelocityVsRelativePhase(self, dTheta, dTheta_t, nBins=100):

        dTheta = dTheta[:-1]

        bins = np.linspace(np.min(dTheta), np.max(dTheta), nBins+1)

        dTheta_t_accum = []

        for i in range(0,nBins):
            dTheta_t_accum.append([])

        for i in range(0, nBins):
            set1 = np.zeros_like(dTheta)
            set2 = np.zeros_like(dTheta)
            set1[np.where(dTheta>bins[i])] = 1
            set2[np.where(dTheta<=bins[i+1])] = 1
            dTheta_t_accum[i] = dTheta_t[np.where(set1*set2==1)[0]]

        dTheta_t_dTheta = np.zeros(len(dTheta_t_accum))
        for i in range(0, len(dTheta_t_accum)):
            if len(dTheta_t_accum[i]) == 0:
                dTheta_t_dTheta[i] = np.nan
            else:
                dTheta_t_dTheta[i] = np.mean(dTheta_t_accum[i])

        dTheta_t_bins = bins[0:-1]

        dTheta_t_bins = np.delete(dTheta_t_bins, np.where(np.isnan(dTheta_t_dTheta)))
        dTheta_t_dTheta = np.delete(dTheta_t_dTheta, np.where(np.isnan(dTheta_t_dTheta)))

        return dTheta_t_bins, dTheta_t_dTheta
                         
            
    #######################################################
        
    def plotRelativePhaseVelocityVsRelativePhase(self):

        self.dTheta_t_bins, self.dTheta_t_dTheta = self.genRelativePhaseVelocityVsRelativePhase(self.dTheta, self.dTheta_t)
        
        plt.plot(self.dTheta_t_bins, self.dTheta_t_dTheta)
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
        
        ax1.plot(self.t2, self.dTheta)#, t2, dtheta2, t3, dtheta3, t4, dtheta4)
    
        ax1.set_ylim([0,1])
        ax1.set_yticks(np.linspace(0,1,3))
        ax1.set_ylabel(r"$\phi_{cell} - \phi_{sub}$", fontsize=16)
        ax1.set_xlabel("t (sec)")
        ax1.set_title(r"$\Delta$=" + '%.2f' % self.Delta, fontsize=16)

        plt.savefig(self.title + 'relativePhase.eps', dpi=160, facecolor='w')
        plt.show()
    





