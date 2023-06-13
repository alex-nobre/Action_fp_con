# -*- coding: utf-8 -*-

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
kList= [i for i in range(1,15)] # max k possible is 15, and need to be integers
rList= np.linspace(-10,-0.1,num=15)#np.linspace(-4,-1,num=10)
cList= np.linspace(0,0.0021,num=15)#np.linspace(0,0.0003,num=10)

startVals=[4., 375.]

# List to store results
avFits = [None] * (len(kList) * len(rList) * len(cList))

fitEl = 0
for ik in range(len(kList)):
  for ir in range(len(rList)):
    for ic in range(len(cList)):
      
      # Generate model
      k=kList[ik]
      r=rList[ir]
      c=cList[ic]
      fmtp = fMTP(r, c, k)
      
      # List for fits by FP
      simList = [None] * len(FPs)
      
      # Run fit by FP 
      for iFP in range(len(FPs)):
        
        FP = FPs[[iFP]]
        
        # Simulate experiment
        exp = FPexp(FPs = FP, distribution = "constant", tr_per_block = 150)
        
        sim, prep = exp.run_exp(fmtp)
        sim = sim.iloc[1:, :] # first trial has no preparation
        
        # AVerage preparation at FP
        sim = sim.groupby(["FP"])[["prep"]].mean().reset_index()
        
        simList[iFP] = sim
      
      
      # Create data frame from sims
      simAv = pd.DataFrame(np.array(simList).reshape(2,2),  columns = ["FP", "prep"])
      
      # Fit data for action and external data 
      sseAction = op.minimize(sse, startVals, args = (dataAction, simAv), method = 'L-BFGS-B')
      sseExternal = op.minimize(sse, startVals, args = (dataExternal, simAv), method = 'L-BFGS-B')
      
      # Get R2 values for each condition
      RsseAction = np.corrcoef((sseAction.x[0] * simAv.prep + sseAction.x[1]), dataAction.RT)[0][1]**2
      RsseExternal = np.corrcoef((sseExternal.x[0] * simAv.prep + sseExternal.x[1]), dataExternal.RT)[0][1]**2
      
      # Store fits for action and external conditions
      # provList[0] = (k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
      # sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
      
      avFits[fitEl] = (k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
      sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
      fitEl += 1

avFitResults = pd.DataFrame(avFits, 
columns = ['k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction', 'R2_SSE_Action', 
'a SSE External', 'b SSE External', 'SSEExternal', 'R2_SSE_External'])

avFitResults.to_csv("./Modeling/fMTP/Fitting/" + "avFitResults.csv")
