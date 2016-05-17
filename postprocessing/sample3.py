import numpy as np
import matplotlib.pyplot as plt
import experimentLib

# Data extracted from stretch sample 3, experiment conducted 20160516

experimentTitle = '20160516_substrate3'

title = ['sample3_21min_0v', 'sample3_24min_0v', 'sample3_27min_0v', 'sample3_30min_1_3v', 'sample3_32min_1_5v', 'sample3_34min_1_7v', 'sample3_36min_1_7v', 'sample3_38min_1_9v', 'sample3_40min_1_9v', 'sample3_43min_2_1v', 'sample3_46min_2_1v', 'sample3_49min_0v']

dt = 1/200.0


#######################################
# Optional:  Use nominal substrate frequency (no phase data!)

useAvailableSubstrateEvents = 1  # Set to unity to use *data* for substrate, as opposed to measured frequency

nCycles = 5
nominalSubPeriod = dt / nCycles * np.array([0, 0, 0, (2228-182), (2061-216), (1586-95), (1561-79), (1429-70), (1510-153), (1300-119), (1297-122), 0])

nominalSubFreq = np.zeros_like(nominalSubPeriod)

for i in range(0,nominalSubPeriod.shape[0]):
    if nominalSubPeriod[i] == 0:
        nominalSubFreq[i] = 0
    else:
        nominalSubFreq[i] = 1/nominalSubPeriod[i]


#################################################
############n#####################################

# Contractions (reaction time frame @ 25fps playback)

reactionTime = dt * 5

c3_21min_0v = dt * np.array([68, 493, 922, 1350, 1782, 2228, 2673, 3112, 3541])#, 3781])
c3_24min_0v = dt * np.array([440, 947, 1437, 1949, 2433, 2900, 3402, 3898])
c3_27min_0v = dt * np.array([19, 497, 971, 1431, 1913, 2403, 2868, 3326, 3810])
c3_30min_1_3v = dt * np.array([74, 499, 909, 1304, 1636, 2026, 2522, 2850, 3263, 3776])
c3_32min_1_5v = dt * np.array([134, 502, 869, 1243, 1610, 1979, 2345, 2714, 3084, 3453, 3818])
c3_34min_1_7v = dt * np.array([224, 793, 1119, 1495, 1981, 2317, 2631, 3017, 3465, 3788, 4100])
c3_36min_1_7v = dt * np.array([17, 312, 615, 907, 1205, 1548, 2095, 2467, 2978, 3363, 3879])
c3_38min_1_9v = dt * np.array([146, 574, 984, 1473, 1811, 2243, 2617, 3161, 3594, 3979])
c3_40min_1_9v = dt * np.array([270, 619, 1059, 1435, 1875, 2367, 2710, 2940, 3156, 3407, 3730, 4043, 4380, 4855, 5271, 5686])
c3_43min_2_1v = dt * np.array([196, 654, 998, 1381, 1740, 2057, 2337, 2729, 3042, 3491, 3829, 4209, 4541, 4919, 5210, 5613, 6108])
c3_46min_2_1v = dt * np.array([450, 702, 1140, 1626, 2087, 2457, 2817, 3262, 3628, 3986, 4411, 4771, 5160, 5556, 5811, 6064])
c3_49min_0v = dt * np.array([123, 507, 918, 1349, 1686, 2104, 2449, 2900])


#################################################
#################################################

# Reference substrate stretch (approximated with unestablished phase reference by max displacement from rigid body stack reg.  Data points are *inconsistent from video to video* !!!):

# subEventTimeShift is the timeshift to set to a desired event. Set such that the *SUM* of timeshift and events yields the time of minimum substrate strain (which will later be mapped to phase value of zero)


s3_21min_0v = dt * np.array([])

s3_24min_0v = dt * np.array([])

s3_27min_0v = dt * np.array([])

s3_30min_1_3v = dt * np.array([])

s3_32min_1_5v = dt * np.array([])

s3_34min_1_7v = dt * np.array([])

s3_36min_1_7v = dt * np.array([80, 377, 674, 970, 1266, 1560, 1857, 2156, 2449, 2748, 3046, 3343, 3638, 3937]) + 52 * dt

s3_38min_1_9v = dt * np.array([])

s3_40min_1_9v = dt * np.array([23, 294, 565, 838, 1108, 1379, 1654, 1923, 2194, 2468, 2738, 3010, 3284, 3557, 3829, 4100, 4371, 4644, 4917, 5187, 5458, 5731, 6003]) + 49*dt

s3_43min_2_1v = dt * np.array([])

s3_46min_2_1v = dt * np.array([])

s3_49min_0v = dt * np.array([])



#################################################
#################################################

voltage = np.array([0, 0, 0, 0, 1.3, 1.5, 1.7, 1.7, 1.9, 1.9, 2.1, 2.1, 0])

startTime = np.array([21, 24, 27, 30, 32, 34, 36, 38, 40, 43, 46, 49])    # In minutes!

cellEvents = [c3_21min_0v, c3_24min_0v, c3_27min_0v, c3_30min_1_3v, c3_32min_1_5v, c3_34min_1_7v, c3_36min_1_7v, c3_38min_1_9v, c3_40min_1_9v, c3_43min_2_1v, c3_46min_2_1v, c3_49min_0v]

subEvents = [s3_21min_0v, s3_24min_0v, s3_27min_0v, s3_30min_1_3v, s3_32min_1_5v, s3_34min_1_7v, s3_36min_1_7v, s3_38min_1_9v, s3_40min_1_9v, s3_43min_2_1v, s3_46min_2_1v, s3_49min_0v]

subEventTimeShift = np.zeros_like(voltage)

#################################################
#################################################

params = {'cellEvents':cellEvents, 'subEvents':subEvents, 'voltage':voltage, 'startTime':startTime, 'dt':dt, 'useAvailableSubstrateEvents':useAvailableSubstrateEvents, 'nominalSubFreq':nominalSubFreq, 'reactionTime':reactionTime, 'subEventTimeShift':subEventTimeShift, 'title':title, 'experimentTitle':experimentTitle}

s3 = experimentLib.experiment(params)

m8 = s3.genMeasurement(8)
m8.hist_substratePhasesOnContraction()


m6 = s3.genMeasurement(6)
m6.hist_substratePhasesOnContraction()
