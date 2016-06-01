import numpy as np
import pdb
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import measurementLib2 as measurementLib
from matplotlib.ticker import FormatStrFormatter




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
        self.maxStrain = params['maxStrain']
        self.cellNaturalFreq = params['cellNaturalFreq']
        self.title = params['title']
        self.experimentTitle = params['experimentTitle']
        
        

        # If nonzero reaction time provided, substract from cellEvents:        
        self.n = len(self.title)
        self.subtractReactionTime()
        self.shiftSubstrateEvents()

        self.computeCellFreq()
        self.computeSubFreq()
        self.computeDelta()

        
    #############################################

    def genStretchedMeasurement(self, i):

        return measurementLib.stretchedMeasurement(self.cellEvents[i], self.cellFreq[i], self.cellNaturalFreq, self.dt, self.startTime[i], self.title[i], self.subEvents[i], self.subFreq[i], self.Delta[i])


    #############################################
    
    def genUnstretchedMeasurement(self, i):
        return measurementLib.unstretchedMeasurement(self.cellEvents[i], self.cellFreq[i],  self.cellNaturalFreq, self.dt, self.startTime[i], self.title[i])

    
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


    #############################################

    def computeDelta(self):
        self.Delta = np.zeros(self.n)

        for i in range(0, self.n):
            if self.subFreq[i] !=0:
                self.Delta[i] = (self.subFreq[i] - self.cellNaturalFreq)/self.cellNaturalFreq
        

    ##############################################

    def subtractReactionTime(self):
        for i in range(0,self.n):
            self.cellEvents[i] = self.cellEvents[i] - self.reactionTime


    ##############################################
            
    def shiftSubstrateEvents(self):
        for i in range(0, self.n):
            self.subEvents[i] = self.subEvents[i] + self.subEventTimeShift[i]
            
        
    ##############################################
    
    def plotFrequencyErrorbars(self, figsize=(4.5,4.5), top=0.9, bottom=0.15, left=0.14, right=0.9, hspace=0.2, wspace=0.3):
        
        fig = plt.figure(figsize=figsize)
        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right)
        ax1 = fig.add_subplot(111)
        ax1.plot(self.startTime, self.cellFreq, 'r.-', label='cell')
        ax1.errorbar(self.startTime, self.cellFreq, yerr=2*self.cellStd, fmt='o')
        ax1.plot(self.startTime, self.subFreq, 'b.-', label='substrate')

        # y labels:
        ax1.set_yticks(np.linspace(0,np.max([np.max(self.cellFreq),np.max(self.subFreq)]),3))
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        
        # X labels:
        
        #labels = self.title
        #ax1.set_xticklabels(labels, rotation='vertical')

        ax1.set_xlim([np.min(self.startTime)-2,np.max(self.startTime) + 2])
        ax1.set_xticks(np.linspace(self.startTime[0],self.startTime[-1],4))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%2d'))
        ax1.set_xlabel('time (min)')


        ax1.set_ylabel('Frequency (Hz)')
        ax1.set_ylim([0, 1.3 * np.max(np.max(self.subFreq), np.max(self.cellFreq))])
        ax1.set_xlim([np.min(self.startTime) - 2, np.max(self.startTime) + 2])
        legend = ax1.legend(loc='upper left', shadow=True)
        plt.savefig(self.experimentTitle + 'Frequencies.eps', dpi=160, facecolor='w')
        plt.show()
        

    ###############################################
    
    def plotHistograms(self, measurementList, nRows, nColumns, figsize, top=0.9, bottom=0.16, left=0.15, right=0.9, hspace=0.2, wspace=0.3,  hist_maxProbability=0.5):

        nMeasurements = len(measurementList)
        fig, axs = plt.subplots(nRows, nColumns, figsize=figsize, facecolor='w', edgecolor='k')
        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace) 
        axs = axs.ravel()
        
        tt = np.linspace(0,2*np.pi,1000)
        yy = self.maxStrain * (0.5 - 0.5 * np.cos(tt))
        nBins = 8
        bins = np.linspace(0,2*np.pi,nBins+1)
        
        for i in range(nMeasurements):
        
            hist, rbins = np.histogram(2*np.pi*measurementList[i].subTheta[measurementList[i].cellIx], bins=bins)
            widths = np.diff(bins)
            scaling = 1/measurementList[i].cellIx.size
            hist = hist*scaling
            axs[i].bar(rbins[:-1], hist, widths)

            if self.cellNaturalFreq == 0:
                print("Don't forget to update cellNaturalFreq")
            axs[i].set_xlim([0,2*np.pi])
            axs[i].set_xticks(np.linspace(0,2*np.pi,5))
            axs[i].set_xticklabels(['0','',r"$\pi$",'',r"2$\pi$"])
            axs[i].set_ylim([0,hist_maxProbability])
            axs[i].set_yticks(np.linspace(0,hist_maxProbability,3))
            axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            #axs[i].set_title(r"$\Delta$=" + str(measurementList[i].Delta), fontsize=16)
            axs[i].set_title(r"$\Delta$=" + '%.2f' % measurementList[i].Delta, fontsize=16)
            
            for tl in axs[i].get_yticklabels():
                tl.set_color('b')
        
            ax2 = axs[i].twinx()
            ax2.plot(tt,yy,'r--')
            ax2.set_yticks(np.linspace(0,self.maxStrain,4))
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            for tl in ax2.get_yticklabels():
                tl.set_color('r')
        
            if np.mod(i,nColumns)==0:
                ax2.set_yticklabels([])
                axs[i].set_ylabel(r"$p_{c}$", color='blue', fontsize=20)
            elif np.mod(i,nColumns)==nColumns-1:
                axs[i].set_yticklabels([])
                ax2.set_ylabel(r"$\varepsilon$", color='red', fontsize=20)
            else:
                axs[i].set_yticklabels([])
                ax2.set_yticklabels([])
        
            if np.floor(i/nColumns)+1==nRows:
                axs[i].set_xlabel(r"$\phi_{substrate}$", fontsize=16)

        plt.savefig(self.experimentTitle + '_histograms.eps', dpi=160, facecolor='w')
        plt.show()

        
        ################################################


    def stretchedTauHistogram(self, measurementList, nRows, nColumns, figsize, top=0.9, bottom=0.16, left=0.15, right=0.9, hspace=0.2, wspace=0.3,  hist_maxProbability=0.5, nBins=16):

        epsilon = 0.02


        nMeasurements = len(measurementList)

        fig, axs = plt.subplots(nRows, nColumns, figsize=figsize, facecolor='w', edgecolor='k')

        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace) 
        axs = axs.ravel()

        tauMax = 0
        for i in range(nMeasurements):
            if np.max(np.diff(measurementList[i].cellEvents)) > tauMax:
                tauMax = np.max(np.diff(measurementList[i].cellEvents))
            elif 1/measurementList[i].subFreq > tauMax:
                tauMax = 1/measurementList[i].subFreq
                
        bins = np.linspace(0,tauMax,nBins+1)
        
        for i in range(nMeasurements):

            axs[i].add_patch(
                patches.Rectangle(
                    (1/measurementList[i].subFreq - epsilon, 0),    # (x,y)
                    2*epsilon,                 # width
                    hist_maxProbability,       # height
                    facecolor="cyan"))

            axs[i].add_patch(
                patches.Rectangle(
                    (1/measurementList[i].cellNaturalFreq - epsilon, 0),
                    2*epsilon,
                    hist_maxProbability,
                    facecolor="red"))
            
            hist, rbins = np.histogram(np.diff(measurementList[i].cellEvents), bins=bins)
            widths = np.diff(bins)
            scaling = 1/(np.diff(measurementList[i].cellEvents).size)
            hist = hist*scaling
            axs[i].bar(rbins[:-1], hist, widths)

            if self.cellNaturalFreq == 0:
                print("Don't forget to update cellNaturalFreq")
                    
            axs[i].set_xlim([0,tauMax])
            axs[i].set_xticks(np.linspace(0,tauMax,5))
            axs[i].set_xticklabels(['0','',tauMax/2,'',tauMax])
            axs[i].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            
            axs[i].set_ylim([0,hist_maxProbability])
            axs[i].set_yticks(np.linspace(0,np.min([hist_maxProbability,1]),3))
            axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            #axs[i].set_title(r"$\Delta$=" + str(measurementList[i].Delta), fontsize=16)
            axs[i].set_title(r"$\Delta$=" + '%.2f' % measurementList[i].Delta, fontsize=16)
            
            
            for tl in axs[i].get_yticklabels():
                tl.set_color('b')
        
        
            if np.mod(i,nColumns)==0:
                axs[i].set_ylabel(r"$p_{c}$", color='blue', fontsize=20)
            elif np.mod(i,nColumns)==nColumns-1:
                axs[i].set_yticklabels([])
            else:
                axs[i].set_yticklabels([])
        
            if np.floor(i/nColumns)+1==nRows:
                axs[i].set_xlabel(r"$\tau$ (s)", fontsize=16)
                
            #axs[i].hist( measurementList[i].subTheta[measurementList[i].cellIx],  bins=bins, normed=True)

        plt.savefig(self.experimentTitle + '_stretchedTauHistograms.eps', dpi=160, facecolor='w')
        plt.show()


    ##########################################################
        

    def unstretchedTauHistogram(self, measurementList, nRows, nColumns, figsize, top=0.9, bottom=0.16, left=0.15, right=0.9, hspace=0.2, wspace=0.3,  hist_maxProbability=0.5, nBins=16):

        epsilon = 0.02

        nMeasurements = len(measurementList)
        fig, axs = plt.subplots(nRows, nColumns, figsize=figsize, facecolor='w', edgecolor='k')
        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace) 
        axs = axs.ravel()

        tauMax = 0
        for i in range(nMeasurements):
            if 1.1*np.max(np.diff(measurementList[i].cellEvents)) > tauMax:
                tauMax = 1.1*np.max(np.diff(measurementList[i].cellEvents))
                
        bins = np.linspace(0,tauMax,nBins+1)
        
        for i in range(nMeasurements):

            axs[i].add_patch(
                patches.Rectangle(
                    (1/measurementList[i].cellNaturalFreq - epsilon, 0),
                    2*epsilon,
                    hist_maxProbability,
                    facecolor="red"))
            
            hist, rbins = np.histogram(np.diff(measurementList[i].cellEvents), bins=bins)
            widths = np.diff(bins)
            scaling = 1/(np.diff(measurementList[i].cellEvents).size)
            hist = hist*scaling
            axs[i].bar(rbins[:-1], hist, widths)

            if self.cellNaturalFreq == 0:
                print("Don't forget to update cellNaturalFreq")
                    
            axs[i].set_xlim([0,tauMax])
            axs[i].set_xticks(np.linspace(0,tauMax,5))
            axs[i].set_xticklabels(['0','',tauMax/2,'',tauMax])
            axs[i].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            
            axs[i].set_ylim([0,hist_maxProbability])
            axs[i].set_yticks(np.linspace(0,np.min([hist_maxProbability,1]),3))
            axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            #axs[i].set_title(r"$\Delta$=" + str(measurementList[i].Delta), fontsize=16)
            axs[i].set_title("Unstretched, t=%2d min" % measurementList[i].startTime, fontsize=16)
            
            
            for tl in axs[i].get_yticklabels():
                tl.set_color('b')
        
        
            if np.mod(i,nColumns)==0:
                axs[i].set_ylabel(r"$p_{c}$", color='blue', fontsize=20)
            elif np.mod(i,nColumns)==nColumns-1:
                axs[i].set_yticklabels([])
            else:
                axs[i].set_yticklabels([])
        
            if np.floor(i/nColumns)+1==nRows:
                axs[i].set_xlabel(r"$\tau$ (s)", fontsize=16)
                
            #axs[i].hist( measurementList[i].subTheta[measurementList[i].cellIx],  bins=bins, normed=True)

        plt.savefig(self.experimentTitle + '_unstretchedTauHistograms.eps', dpi=160, facecolor='w')
        plt.show()


    ##########################################################
        

    def plot_dTheta_t_of_dTheta(self, measurementList, min_dTheta_t, max_dTheta_t, nRows, nColumns, figsize, top=0.9, bottom=0.16, left=0.15, right=0.9, hspace=0.2, wspace=0.3):

        nMeasurements = len(measurementList)

        fig, axs = plt.subplots(nRows, nColumns, figsize=figsize, facecolor='w', edgecolor='k')
        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace) 
        axs = axs.ravel()

        for i in range(nMeasurements):

            axs[i].plot(measurementList[i].dTheta_t_bins, measurementList[i].dTheta_t_dTheta)

            axs[i].set_xlim([np.min(measurementList[i].dTheta_t_bins), np.max(measurementList[i].dTheta_t_bins)])

            #axs[i].set_xticks(np.linspace(0,tauMax,5))
            #axs[i].set_xticklabels(['0','',tauMax/2,'',tauMax])
            #axs[i].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            
            axs[i].set_ylim([min_dTheta_t, max_dTheta_t])
            #axs[i].set_yticks(np.linspace(0,np.min([hist_maxProbability,1]),3))
            #axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            #axs[i].set_title(r"$\Delta$=" + str(measurementList[i].Delta), fontsize=16)
            #axs[i].set_title("Unstretched, t=%2d min" % measurementList[i].startTime, fontsize=16)
            
            
            #
            #if np.mod(i,nColumns)==0:
            #    axs[i].set_ylabel(r"$p_{c}$", color='blue', fontsize=20)
            #elif np.mod(i,nColumns)==nColumns-1:
            #    axs[i].set_yticklabels([])
            #else:
            #    axs[i].set_yticklabels([])
            #
            #if np.floor(i/nColumns)+1==nRows:
            #    axs[i].set_xlabel(r"$\tau$ (s)", fontsize=16)
            #    
            #axs[i].hist( measurementList[i].subTheta[measurementList[i].cellIx],  bins=bins, normed=True)

        plt.savefig(self.experimentTitle + '_dTheta_t_of_dTheta.eps', dpi=160, facecolor='w')
        plt.show()







