import numpy as np
import pdb
import matplotlib.pyplot as plt
import measurementLib




class experiment(object):

    def __init__(self, params):

        # Recover parameters:
        self.cellEvents = params['cellEvents']
        self.subEvents = params['subEvents'] 
        self.voltage = params['voltage']
        self.startTime = params['startTime']
        self.dt = params['dt']
        self.useAvailableSubstrateEvents = params['useAvailableSubstrateEvents']
        self.nominalSubFreq = params['nominalSubFreq']
        self.reactionTime = params['reactionTime']
        self.subEventTimeShift = params['subEventTimeShift']
        self.title = params['title']
        self.experimentTitle = params['experimentTitle']
        

        # If nonzero reaction time provided, substract from cellEvents:        
        self.n = len(self.title)
        self.subtractReactionTime()
        self.shiftSubstrateEvents()

        self.computeCellFreq()
        self.computeSubFreq()

        
    #############################################

    def genMeasurement(self, i):
        return measurementLib.measurement(self.cellEvents[i], self.subEvents[i], self.cellFreq[i], self.subFreq[i], self.dt)

    
    #############################################
    
    def computeCellFreq(self):
        self.cellPeriods = []
        self.cellFreq = np.zeros(self.n)
        self.cellStd = np.zeros_like(self.cellFreq)

        for i in range(0, self.n):
            self.cellPeriods.append(np.diff(self.cellEvents[i]))
            self.cellFreq[i] = 1/np.mean(self.cellPeriods[i])
            self.cellStd[i] = np.std(1/self.cellPeriods[i])
        

    #############################################
            
    def computeSubFreq(self):
        self.subFreq = self.nominalSubFreq

        if self.useAvailableSubstrateEvents == 1:
            self.subPeriods = []

            for i in range(0, self.n):
                self.subPeriods.append(np.diff(self.subEvents[i]))
                if self.subPeriods[i].size != 0:
                    self.subFreq[i] = 1/np.mean(self.subPeriods[i])
            

    ##############################################

    def subtractReactionTime(self):
        for i in range(0,self.n):
            self.cellEvents[i] = self.cellEvents[i] - self.reactionTime


    ##############################################
            
    def shiftSubstrateEvents(self):
        for i in range(0, self.n):
            self.subEvents[i] = self.subEvents[i] + self.subEventTimeShift[i]
            
        
    ##############################################
    
    def plotFrequencyErrorbars(self):

        fig = plt.figure(figsize=(6,6))
        fig.subplots_adjust(top=0.98, bottom=0.18, left = 0.14, right = 0.95)
        ax1 = fig.add_subplot(111)
        ax1.plot(self.startTime, self.cellFreq, 'r.-', label='cell')
        ax1.errorbar(self.startTime, self.cellFreq, yerr=2*self.cellStd, fmt='o')
        ax1.plot(self.startTime, self.subFreq, 'b.-', label='substrate')
        ax1.set_xticks(self.startTime)
        labels = self.title
        ax1.set_xticklabels(labels, rotation='vertical')
        ax1.set_ylabel('Frequency (Hz)')
        ax1.set_ylim([0, 1.5 * np.max(np.max(self.subFreq), np.max(self.cellFreq))])
        ax1.set_xlim([np.min(self.startTime) - 2, np.max(self.startTime) + 2])
        legend = ax1.legend(loc='upper right', shadow=True)
        plt.savefig(self.experimentTitle + '.eps', dpi=160, facecolor='w')
        plt.show()
        
        







