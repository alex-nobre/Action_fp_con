# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns

# Import classes
import sys
sys.path.insert(0, "./Modeling/fMTP/Fitting") # path to where classes are stored

from fmtp import fMTP, FPexp

#========================= Functions =======================
# Function to transform preparation values to RT
def temp2RT (prepValue, a, b):
  return (a * prepValue + b)

#====================== Analyze fit ==============================
# Read empirical data and dataframe with fit results
avFitResults = pd.read_csv("./Modeling/fMTP/Fitting/avFitResults.csv")
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT":"RT"})

# Summarise RT by FP and condition and discard other variables
empData = empData.groupby(['foreperiod', 'condition'], as_index=False)[['RT']].mean()

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index(drop=True)
dataExternal = empData.loc[empData.condition == "external"].reset_index(drop=True)

# Constants for simulations in both conditions
FPs = np.array([1.0, 2.8])
startVals=[4., 375.]

# Find minimum SSE for each condition
coefsAction = avFitResults[avFitResults.SSEAction == avFitResults.SSEAction.min()]
coefsAction = coefsAction.iloc[[0]]
coefsExternal = avFitResults[avFitResults.SSEExternal == avFitResults.SSEExternal.min()]

# Parameters for each condition
kAction = coefsAction.iloc[0]["k"]
rAction = coefsAction.iloc[0]["r"]
cAction = coefsAction.iloc[0]["c"]
fmtpAction = fMTP(rAction, cAction, kAction)

kExternal = coefsExternal.iloc[0]["k"]
rExternal = coefsExternal.iloc[0]["r"]
cExternal = coefsExternal.iloc[0]["c"]
fmtpExternal = fMTP(rExternal, cExternal, kExternal)

# List to store values by FP
simActionList = [None] * len(FPs)
simExternalList = [None] * len(FPs)
      
# Run fit by FP 
for iFP in range(len(FPs)):
        
  FP = FPs[[iFP]]
        
  # Simulate experiment
  exp = FPexp(FPs = FP, distribution = "constant", tr_per_block = 150)
  
  # Predictions for action
  state_discr_action, state_con_action = exp.run_exp(fmtpAction)
        
  state_discr_action=state_discr_action[1:] # first trial has no preparation
  mean_state_action=state_discr_action.groupby(['FP'])[["prep"]].mean().reset_index()
  mean_state_action['RT'] = temp2RT(mean_state_action.iloc[0]["prep"], coefsAction.iloc[0]["a SSE Action"], coefsAction.iloc[0]["b SSE Action"])
  
  # Predictions for external
  state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
        
  state_discr_external=state_discr_external[1:] # first trial has no preparation
  mean_state_external=state_discr_external.groupby(['FP'])[["prep"]].mean().reset_index()
  mean_state_external['RT'] = temp2RT(mean_state_external.iloc[0]["prep"], coefsExternal.iloc[0]["a SSE External"], coefsExternal.iloc[0]["b SSE External"])
        
  # AVerage preparation at FP
  simActionList[iFP] = mean_state_action
  simExternalList[iFP] = mean_state_external

# Transform to data frame  
simAction = pd.DataFrame(np.array(simActionList).reshape(2,3),  columns = ["FP", "prep", "RT"])

simAction = simAction[["FP", "RT"]]
simAction.FP = (simAction.FP * 1000)
simAction.FP = simAction.FP.astype(int) # Change to int to match data


simExternal = pd.DataFrame(np.array(simExternalList).reshape(2,3),  columns = ["FP", "prep", "RT"])

simExternal = simExternal[["FP", "RT"]]
simExternal.FP = (simExternal.FP * 1000)
simExternal.FP = simExternal.FP.astype(int) # Change to int to match data


# Add columns to merge sim and data
# simAction["condition"] =["action", "action"]
# simAction["valueType"] = ["fMTP", "fMTP"]
# 
# dataAction["valueType"] = ["data", "data"]
# dataAction = dataAction[["foreperiod", "RT", "condition", "valueType"]] # Reorder to match sim DF
# 
# actionDF = pd.DataFrame(pd.concat([simAction, dataAction.rename(columns = {"foreperiod" : "FP"})]))


# Add columns to merge sim and data
# simExternal["condition"] =["external", "external"]
# simExternal["valueType"] = ["fMTP", "fMTP"]
# 
# dataExternal["valueType"] = ["data", "data"]
# dataExternal = dataExternal[["foreperiod", "RT", "condition", "valueType"]] # Reorder to match sim DF
# 
# externalDF = pd.DataFrame(pd.concat([simExternal, dataExternal.rename(columns = {"foreperiod" : "FP"})]))
# 
# fitDF = pd.DataFrame(pd.concat([actionDF, externalDF])).reset_index(drop = True)
# dataDF = fitDF.loc[fitDF.valueType == "data"].reset_index(drop = True)
# simDF = fitDF.loc[fitDF.valueType == "fMTP"].reset_index(drop = True)

# Rename foreperiod columns in data DFs
dataAction = dataAction.rename(columns = {"foreperiod":"FP"})
dataExternal = dataExternal.rename(columns = {"foreperiod":"FP"})

# Plot
fig1 = plt.figure()
ax = fig1.add_subplot(1,1,1)

#sns.pointplot(x = "FP", y = "RT", hue = "condition")
#sns.catplot(x = "FP", y = "RT", kind = "point", hue = "condition", data = fitDF, col = "valueType")
ax.plot(simAction.FP, simAction.RT, ".--", color = "blue")
ax.plot(dataAction.FP, dataAction.RT, ".-", color = "blue")
ax.plot(simExternal.FP, simExternal.RT, ".--", color = "orange")
ax.plot(dataExternal.FP, dataExternal.RT, ".-", color = "orange")


plt.show()
