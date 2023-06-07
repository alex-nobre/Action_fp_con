
# Import libraries
import numpy as np
import pandas as pd
import os
from scipy import optimize as op
from scipy import stats

# Plotting
import matplotlib.pyplot as plt

# Import fmtp classes
import sys
sys.path.insert(0, './Modeling/fMTP/Fitting') # path to where classes are stored

from fmtp import fMTP, FPexp
from hazard import fMTPhz
from fit import sort_fit, show_fit, get_fit 

#=================================== functions =============================

def temp2RT (prepValue, a, b):
    return (a*prepValue + b)

def sse (parm, emp, sim):
    simRTs = temp2RT(sim.prep, parm[0], parm[1])
    SSE = sum(pow(emp.RT - simRTs, 2))
    return SSE

def rmse (parm, emp, sim):
    simRTs = temp2RT(sim.prep, parm[0], parm[1])
    simzRTs = stats.zscore(simRTs)
    RMSE = np.sqrt(((emp.zRT - simzRTs)**2).mean())
    return RMSE

def anticorr (parm, emp, sim):
    simRTs = temp2RT(sim.prep, parm[0], parm[1])
    antiCorr = 1 - np.corrcoef(emp.RT, simRTs)[0][1]
    return antiCorr

#============================== Fit hazard =============================    

# Read empirical data
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT": "RT"})

# FP values
FPs = np.array([1.0, 2.8])

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index()
dataExternal = empData.loc[empData.condition == "external"].reset_index()

# Summarise RT by FP and discard other variables
dataAction = dataAction.groupby(["foreperiod"], as_index = False)[["RT"]].mean()
dataExternal = dataExternal.groupby(["foreperiod"], as_index = False)[["RT"]].mean()


# Set paramters for fmtp
k = 4 # temporal precision
r = -2.81 # rate of forgettting
c = 0.0001 # memory persistence
fmtp = fMTP(r, c, k)

#FP = FPs[0]
startVals=[4., 375.]

simList = [None] * len(FPs)

for iFP in range(len(FPs)):
  
  FP = FPs[[iFP]]
  
  # Run fit by FP 
  exp = FPexp(FPs = x, distribution = "constant", tr_per_block = 150)
  
  sim, prep = exp.run_exp(fmtp)
  sim = sim.iloc[1:, :] # first trial has no preparation
  
  # AVerage preparation at FP
  sim = sim.groupby(["FP"])[["prep"]].mean().reset_index()
  
  simList[iFP] = sim


# Create data frame from sims
simAv = pd.DataFrame(np.array(simList).reshape(2,2),  columns = ["FP", "prep"])


# Compare fits for action and external conditions
sseAction = op.minimize(sse, startVals, args = (dataAction, simAv), method = 'L-BFGS-B')
sseExternal = op.minimize(sse, startVals, args = (dataExternal, simAv), method = 'L-BFGS-B')


# Fit data for action and external data  

