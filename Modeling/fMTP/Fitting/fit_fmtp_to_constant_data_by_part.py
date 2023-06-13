# -*- coding: utf-8 -*-

reticulate::repl_python()

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

#=================================== Functions =============================

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

#============================================================================#
#========================= 1. All free parameters  ===========================    
#============================================================================#

# Read empirical data
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT": "RT"})

IDList=np.unique(empData['ID']).tolist()

# FP values
FPs = np.array([1.0, 2.8])

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index()
dataExternal = empData.loc[empData.condition == "external"].reset_index()

# Summarise RT by FP and discard other variables
dataAction = dataAction.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()
dataExternal = dataExternal.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()


# Set paramters for fmtp
kList= [i for i in range(1,15)] # max k possible is 15, and need to be integers
rList= np.linspace(-10,-0.1,num=15) 
cList= np.linspace(0,0.0021,num=15)

startVals=[4., 375.]

# List to store results
subFits = [None] * (len(kList) * len(rList) * len(cList) *len(IDList))

fitEl = 0
for i, iID in enumerate(IDList):
  print(i)
  
  # Action and external datasets for part
  IDDataAction=dataAction.loc[dataAction.ID==iID].reset_index(drop=True)
  IDDataExternal=dataExternal.loc[dataExternal.ID==iID].reset_index(drop=True)
  
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
        simDF = pd.DataFrame(np.array(simList).reshape(2,2),  columns = ["FP", "prep"])
        
        # Fit data for action and external data 
        sseAction = op.minimize(sse, startVals, args = (IDDataAction, simDF), method = 'L-BFGS-B')
        sseExternal = op.minimize(sse, startVals, args = (IDDataExternal, simDF), method = 'L-BFGS-B')
        
        # Get R2 values for each condition
        RsseAction = np.corrcoef((sseAction.x[0] * simDF.prep + sseAction.x[1]), IDDataAction.RT)[0][1]**2
        RsseExternal = np.corrcoef((sseExternal.x[0] * simDF.prep + sseExternal.x[1]), IDDataExternal.RT)[0][1]**2
        
        # Store fits for action and external conditions
        # provList[0] = (k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
        # sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
        
        subFits[fitEl] = (iID, k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
        sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
        fitEl += 1

subFitResults = pd.DataFrame(subFits, 
columns = ["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action", 
"a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"])

subFitResults.to_csv("./Modeling/fMTP/Fitting/" + "subFitResults.csv")

#===============================================================================#
#============================ 2. c and r fixed ==================================
#===============================================================================#
# Read empirical data
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT": "RT"})

IDList=np.unique(empData['ID']).tolist()

# FP values
FPs = np.array([1.0, 2.8])

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index()
dataExternal = empData.loc[empData.condition == "external"].reset_index()

# Summarise RT by FP and discard other variables
dataAction = dataAction.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()
dataExternal = dataExternal.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()

# Set paramters for fmtp
kList= [i for i in range(1,15)] # max k possible is 15, and need to be integers
r = -2.81 # rate of forgettting
c = 0.0001 # memory persistence

startVals=[4., 375.]

# List to store results
subFits = [None] * (len(kList) *len(IDList))

fitEl = 0
for i, iID in enumerate(IDList):
  print(i)
  
  # Action and external datasets for part
  IDDataAction=dataAction.loc[dataAction.ID==iID].reset_index(drop=True)
  IDDataExternal=dataExternal.loc[dataExternal.ID==iID].reset_index(drop=True)
  
  for ik in range(len(kList)):
        
    # Generate model
    k=kList[ik]
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
    simDF = pd.DataFrame(np.array(simList).reshape(2,2),  columns = ["FP", "prep"])
    
    # Fit data for action and external data 
    sseAction = op.minimize(sse, startVals, args = (IDDataAction, simDF), method = 'L-BFGS-B')
    sseExternal = op.minimize(sse, startVals, args = (IDDataExternal, simDF), method = 'L-BFGS-B')
    
    # Get R2 values for each condition
    RsseAction = np.corrcoef((sseAction.x[0] * simDF.prep + sseAction.x[1]), IDDataAction.RT)[0][1]**2
    RsseExternal = np.corrcoef((sseExternal.x[0] * simDF.prep + sseExternal.x[1]), IDDataExternal.RT)[0][1]**2
    
    # Store fits for action and external conditions
    # provList[0] = (k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
    # sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
    
    subFits[fitEl] = (iID, k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
    sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
    fitEl += 1

subFitResults = pd.DataFrame(subFits, 
columns = ["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action", 
"a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"])

subFitResults.to_csv("./Modeling/fMTP/Fitting/" + "subFitResults_k.csv")


#===============================================================================#
#============================ 3. k and c fixed ==================================
#===============================================================================#
# Read empirical data
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT": "RT"})

IDList=np.unique(empData['ID']).tolist()

# FP values
FPs = np.array([1.0, 2.8])

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index()
dataExternal = empData.loc[empData.condition == "external"].reset_index()

# Summarise RT by FP and discard other variables
dataAction = dataAction.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()
dataExternal = dataExternal.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()

# Set paramters for fmtp
k = 4 # temporal precision
rList= np.linspace(-10,-0.1,num=15) # rate of forgettting
c = 0.0001 # memory persistence

startVals=[4., 375.]

# List to store results
subFits = [None] * (len(rList) *len(IDList))

fitEl = 0
for i, iID in enumerate(IDList):
  print(i)
  
  # Action and external datasets for part
  IDDataAction=dataAction.loc[dataAction.ID==iID].reset_index(drop=True)
  IDDataExternal=dataExternal.loc[dataExternal.ID==iID].reset_index(drop=True)
  
  for ir in range(len(rList)):
        
    # Generate model
    r=rList[ir]
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
    simDF = pd.DataFrame(np.array(simList).reshape(2,2),  columns = ["FP", "prep"])
    
    # Fit data for action and external data 
    sseAction = op.minimize(sse, startVals, args = (IDDataAction, simDF), method = 'L-BFGS-B')
    sseExternal = op.minimize(sse, startVals, args = (IDDataExternal, simDF), method = 'L-BFGS-B')
    
    # Get R2 values for each condition
    RsseAction = np.corrcoef((sseAction.x[0] * simDF.prep + sseAction.x[1]), IDDataAction.RT)[0][1]**2
    RsseExternal = np.corrcoef((sseExternal.x[0] * simDF.prep + sseExternal.x[1]), IDDataExternal.RT)[0][1]**2
    
    # Store fits for action and external conditions
    # provList[0] = (k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
    # sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
    
    subFits[fitEl] = (iID, k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
    sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
    fitEl += 1

subFitResults = pd.DataFrame(subFits, 
columns = ["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action", 
"a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"])

subFitResults.to_csv("./Modeling/fMTP/Fitting/" + "subFitResults_r.csv")


#===============================================================================#
#============================ 4. k and r fixed ==================================
#===============================================================================#
# Read empirical data
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT": "RT"})

IDList=np.unique(empData['ID']).tolist()

# FP values
FPs = np.array([1.0, 2.8])

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index()
dataExternal = empData.loc[empData.condition == "external"].reset_index()

# Summarise RT by FP and discard other variables
dataAction = dataAction.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()
dataExternal = dataExternal.groupby(["ID", "foreperiod"], as_index = False)[["RT"]].mean()

# Set paramters for fmtp
k = 4 # temporal precision
r = -2.81 # rate of forgettting
cList= np.linspace(0,0.0021,num=15) # memory persistence

startVals=[4., 375.]

# List to store results
subFits = [None] * (len(cList) *len(IDList))

fitEl = 0
for i, iID in enumerate(IDList):
  print(i)
  
  # Action and external datasets for part
  IDDataAction=dataAction.loc[dataAction.ID==iID].reset_index(drop=True)
  IDDataExternal=dataExternal.loc[dataExternal.ID==iID].reset_index(drop=True)
  
  for ic in range(len(cList)):
        
    # Generate model
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
    simDF = pd.DataFrame(np.array(simList).reshape(2,2),  columns = ["FP", "prep"])
    
    # Fit data for action and external data 
    sseAction = op.minimize(sse, startVals, args = (IDDataAction, simDF), method = 'L-BFGS-B')
    sseExternal = op.minimize(sse, startVals, args = (IDDataExternal, simDF), method = 'L-BFGS-B')
    
    # Get R2 values for each condition
    RsseAction = np.corrcoef((sseAction.x[0] * simDF.prep + sseAction.x[1]), IDDataAction.RT)[0][1]**2
    RsseExternal = np.corrcoef((sseExternal.x[0] * simDF.prep + sseExternal.x[1]), IDDataExternal.RT)[0][1]**2
    
    # Store fits for action and external conditions
    # provList[0] = (k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
    # sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
    
    subFits[fitEl] = (iID, k, r, c, sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction,
    sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal)
    fitEl += 1

subFitResults = pd.DataFrame(subFits, 
columns = ["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action", 
"a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"])

subFitResults.to_csv("./Modeling/fMTP/Fitting/" + "subFitResults_c.csv")
