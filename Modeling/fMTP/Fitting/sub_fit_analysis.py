# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
from scipy import stats

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns

# Import classes
import sys
sys.path.insert(0, "./Modeling/fMTP/Fitting") # path to where classes are stored

from fmtp import fMTP, FPexp

#=================================== Functions =============================
# Function to transform preparation values to RT
def temp2RT (prepValue, a, b):
  return (a * prepValue + b)

#============================== Analyze fit =============================
# Read empirical data and dataframe with fit results
subFitResults = pd.read_csv("./Modeling/fMTP/Fitting/subFitResults.csv")
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT":"RT"})
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index(drop=True)
dataExternal = empData.loc[empData.condition == "external"].reset_index(drop=True)

# Constants for simulations in both conditions
FPs = np.array([1.0, 2.8])
startVals=[4., 375.]

# Find minimum SSE for each condition
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop=True).drop(columns = ["Unnamed: 0"])
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop=True).drop(columns = ["Unnamed: 0"])

# bestfCoefs = pd.DataFrame(pd.concat([coefsAction.rename(columns = {"ID":"ID", "k":"kAction", "r":"rAction", "c":"cAction", "a SSE Action":"a SSE Action", 
#                                                                    "b SSE Action":"b SSE Action", "SSEAction":"SSE_Action", "R2_SSE_Action":"R2_SSE_Action"}), 
#                                      coefsExternal.rename(columns = {"k":"kExternal", "r":"rExternal", "c":"cExternal", "a SSE External":"a SSE External", 
#                                                                      "b SSE External":"b SSE External", "SSEExternal": "SSE_External", "R2_SSE_External":"R2_SSE_External"}).drop(columns = ["ID"])], 
#                                     axis = 1), columns = ["ID", "kAction", "rAction", "cAction", "a SSE Action", "b SSE Action", "SSE_Action", "R2_SSE_Action", 
#                                                           "kExternal", "rExternal", "cExternal", "a SSE External", "b SSE External", "SSE_External", "R2_SSE_External"])

# Prepare panels for plotting
fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(4,6,sharey=False)
ax = ax.ravel()

titles = coefsAction.ID.tolist()

for iSub, sub in enumerate(titles):
  
  # Parameters for each condition for subject
  subCoefsAction = coefsAction[coefsAction.ID == sub]
  subCoefsExternal = coefsExternal[coefsExternal.ID == sub]
  
  kAction = subCoefsAction.iloc[0]["k"]
  rAction = subCoefsAction.iloc[0]["r"]
  cAction = subCoefsAction.iloc[0]["c"]
  fmtpAction = fMTP(rAction, cAction, kAction)
  
  kExternal = subCoefsExternal.iloc[0]["k"]
  rExternal = subCoefsExternal.iloc[0]["r"]
  cExternal = subCoefsExternal.iloc[0]["c"]
  fmtpExternal = fMTP(rExternal, cExternal, kExternal)
  
  # Subset participant's data by condition
  dataSubAction = dataAction[dataAction['ID']==sub]
  dataSubExternal = dataExternal[dataExternal['ID']==sub]
  
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
    mean_state_action['RT'] = temp2RT(mean_state_action.iloc[0]["prep"], subCoefsAction.iloc[0]["a SSE Action"], subCoefsAction.iloc[0]["b SSE Action"])
    
    # Predictions for external
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
          
    state_discr_external=state_discr_external[1:] # first trial has no preparation
    mean_state_external=state_discr_external.groupby(['FP'])[["prep"]].mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.iloc[0]["prep"], subCoefsExternal.iloc[0]["a SSE External"], subCoefsExternal.iloc[0]["b SSE External"])
          
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
  dataSubAction = dataSubAction.rename(columns = {"foreperiod":"FP"})
  dataSubExternal = dataSubExternal.rename(columns = {"foreperiod":"FP"})
  
  # Plot
  ax[iSub].plot(simAction.FP, simAction.RT, ".--", color = "blue")
  ax[iSub].scatter(dataSubAction.FP, dataSubAction.RT, marker = ".", color = "blue")
  ax[iSub].plot(simExternal.FP, simExternal.RT, ".--", color = "orange")
  ax[iSub].scatter(dataSubExternal.FP, dataSubExternal.RT, marker = ".", color = "orange")
  ax[iSub].set_title(sub)


plt.show()

figname = "./Modeling/fMTP/Fitting/Plots/con_FP_sse_fits.png"
plt.savefig(figname, format = "png", bbox_inches = "tight")


#======================= Compare parameters between conditions ======================
# Build dataset in long format for plotting
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop = True)
coefsAction = coefsAction[["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action"]]
coefsAction["condition"] = pd.Series(["action"] * len(coefsAction))
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop = True)
coefsExternal = coefsExternal[["ID", "k", "r", "c", "a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"]]
coefsExternal["condition"] = pd.Series(["external"] * len(coefsExternal))

coefsLong = pd.DataFrame(pd.concat([coefsAction.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                   "a SSE Action":"aSSE", "b SSE Action":"bSSE", 
                                                                   "SSEAction":"SSE",
                                                                   "R2_SSE_Action":"R2_SSE"}), 
                                   coefsExternal.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                     "a SSE External":"aSSE", 
                                                                     "b SSE External":"bSSE", 
                                                                     "SSEExternal": "SSE",
                                                                     "R2_SSE_External":"R2_SSE"})], 
                                  axis = 0), columns = ["ID", "condition", "k", "r", "c", "aSSE", "bSSE", "SSE", "R2_SSE"]).reset_index(drop = True)


paramFig = plt.figure()
ax = paramFig.subplots(1,3)

sns.stripplot(x = "condition", y = "k", data = coefsLong, jitter = 0.1, orient = "v", ax = ax[0])
sns.stripplot(x = "condition", y = "r", data = coefsLong, jitter = 0.1, orient = "v", ax = ax[1])
sns.stripplot(x = "condition", y = "c", data = coefsLong, jitter = 0.1, orient = "v", ax = ax[2])

plt.show()

paramFig = plt.figure()
ax = paramFig.subplots(3,1)

sns.stripplot(x = "k", y = "condition", data = coefsLong, jitter = 0.1, orient = "h", ax = ax[0])
sns.stripplot(x = "r", y = "condition", data = coefsLong, jitter = 0.1, orient = "v", ax = ax[1])
sns.stripplot(x = "c", y = "condition", data = coefsLong, jitter = 0.1, orient = "v", ax = ax[2])

plt.show()


paramFig = plt.figure()
ax = paramFig.subplots(3,1)

sns.violinplot(x = "k", y = "condition", data = coefsLong, hue = "condition", ax = ax[0])
sns.violinplot(x = "r", y = "condition", data = coefsLong, hue = "condition", ax = ax[1])
sns.violinplot(x = "c", y = "condition", data = coefsLong, hue = "condition", ax = ax[2])

plt.show()

plt.figure()
sns.scatterplot(x = "k", y = "r", data = coefsLong, hue = "condition")
plt.show()

paramFig, ax = plt.subplots(1,1)
colorsList = ["blue"] * len(coefsLong[coefsLong.condition=="action"]) + ["orange"] * len(coefsLong[coefsLong.condition=="external"])
alphaList = ((coefsLong["c"] - min(coefsLong["c"]))/(max(coefsLong["c"]) - min(coefsLong["c"])) * 0.5) + 0.5

# for i in range(len(coefsLong["k"])):
#   ax.scatter(coefsLong["k"][i], coefsLong["r"][i], c = colorsList[i],  s = (coefsLong["c"][i]) * 10000, alpha = alphaList[i],
#   data = coefsLong)
# plt.show()

# Paired samples t-test
stats.ttest_rel(coefsAction.k, coefsExternal.k)
stats.ttest_rel(coefsAction.r, coefsExternal.r)
stats.ttest_rel(coefsAction.c, coefsExternal.c)

#===============================================================================#
#============================ 2. c and r fixed ==================================
#===============================================================================#

# Read empirical data and dataframe with fit results
subFitResults = pd.read_csv("./Modeling/fMTP/Fitting/subFitResults_k.csv")
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT":"RT"})
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index(drop=True)
dataExternal = empData.loc[empData.condition == "external"].reset_index(drop=True)

# Constants for simulations in both conditions
FPs = np.array([1.0, 2.8])
startVals=[4., 375.]

# Find minimum SSE for each condition
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop=True)
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop=True)

# Prepare panels for plotting
fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(4,6,sharey=False)
ax = ax.ravel()

titles = coefsAction.ID.tolist()

for iSub, sub in enumerate(titles):
  
  # Parameters for each condition for subject
  subCoefsAction = coefsAction[coefsAction.ID == sub]
  subCoefsExternal = coefsExternal[coefsExternal.ID == sub]
  
  kAction = subCoefsAction.iloc[0]["k"]
  rAction = subCoefsAction.iloc[0]["r"]
  cAction = subCoefsAction.iloc[0]["c"]
  fmtpAction = fMTP(rAction, cAction, kAction)
  
  kExternal = subCoefsExternal.iloc[0]["k"]
  rExternal = subCoefsExternal.iloc[0]["r"]
  cExternal = subCoefsExternal.iloc[0]["c"]
  fmtpExternal = fMTP(rExternal, cExternal, kExternal)
  
  # Subset participant's data by condition
  dataSubAction = dataAction[dataAction['ID']==sub]
  dataSubExternal = dataExternal[dataExternal['ID']==sub]
  
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
    mean_state_action['RT'] = temp2RT(mean_state_action.iloc[0]["prep"], subCoefsAction.iloc[0]["a SSE Action"], subCoefsAction.iloc[0]["b SSE Action"])
    
    # Predictions for external
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
          
    state_discr_external=state_discr_external[1:] # first trial has no preparation
    mean_state_external=state_discr_external.groupby(['FP'])[["prep"]].mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.iloc[0]["prep"], subCoefsExternal.iloc[0]["a SSE External"], subCoefsExternal.iloc[0]["b SSE External"])
          
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
  
  
  # Rename foreperiod columns in data DFs
  dataSubAction = dataSubAction.rename(columns = {"foreperiod":"FP"})
  dataSubExternal = dataSubExternal.rename(columns = {"foreperiod":"FP"})
  
  # Plot
  ax[iSub].plot(simAction.FP, simAction.RT, ".--", color = "blue")
  ax[iSub].scatter(dataSubAction.FP, dataSubAction.RT, marker = ".", color = "blue")
  ax[iSub].plot(simExternal.FP, simExternal.RT, ".--", color = "orange")
  ax[iSub].scatter(dataSubExternal.FP, dataSubExternal.RT, marker = ".", color = "orange")
  ax[iSub].set_title(sub)


#plt.show()

figname = "./Modeling/fMTP/Fitting/Plots/con_FP_sse_fits_k.png"
plt.savefig(figname, format = "png", bbox_inches = "tight")

#======================= Compare parameters between conditions ======================
# Build dataset in long format for plotting
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop = True)
coefsAction = coefsAction[["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action"]]
coefsAction["condition"] = pd.Series(["action"] * len(coefsAction))
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop = True)
coefsExternal = coefsExternal[["ID", "k", "r", "c", "a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"]]
coefsExternal["condition"] = pd.Series(["external"] * len(coefsExternal))

coefsLong = pd.DataFrame(pd.concat([coefsAction.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                   "a SSE Action":"aSSE", "b SSE Action":"bSSE", 
                                                                   "SSEAction":"SSE",
                                                                   "R2_SSE_Action":"R2_SSE"}), 
                                   coefsExternal.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                     "a SSE External":"aSSE", 
                                                                     "b SSE External":"bSSE", 
                                                                     "SSEExternal": "SSE",
                                                                     "R2_SSE_External":"R2_SSE"})], 
                                  axis = 0), columns = ["ID", "condition", "k", "r", "c", "aSSE", "bSSE", "SSE", "R2_SSE"]).reset_index(drop = True)


# Scatter plot
paramFig, ax = plt.subplots(1,1)

sns.stripplot(x = "condition", y = "k", data = coefsLong, jitter = 0.1, orient = "v")

plt.show()

# Histogram
paramFig, ax = plt.subplots(1,1)

sns.histplot(x = "k", hue = "condition", data = coefsLong)

plt.show()

# Paired samples t-test
stats.ttest_rel(coefsAction.k, coefsExternal.k)


#===============================================================================#
#============================ 3. k and c fixed ==================================
#===============================================================================#

# Read empirical data and dataframe with fit results
subFitResults = pd.read_csv("./Modeling/fMTP/Fitting/subFitResults_r.csv")
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT":"RT"})
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index(drop=True)
dataExternal = empData.loc[empData.condition == "external"].reset_index(drop=True)

# Constants for simulations in both conditions
FPs = np.array([1.0, 2.8])
startVals=[4., 375.]

# Find minimum SSE for each condition
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop=True)
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop=True)

# Prepare panels for plotting
fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(4,6,sharey=False)
ax = ax.ravel()

titles = coefsAction.ID.tolist()

for iSub, sub in enumerate(titles):
  
  # Parameters for each condition for subject
  subCoefsAction = coefsAction[coefsAction.ID == sub]
  subCoefsExternal = coefsExternal[coefsExternal.ID == sub]
  
  kAction = subCoefsAction.iloc[0]["k"]
  rAction = subCoefsAction.iloc[0]["r"]
  cAction = subCoefsAction.iloc[0]["c"]
  fmtpAction = fMTP(rAction, cAction, kAction)
  
  kExternal = subCoefsExternal.iloc[0]["k"]
  rExternal = subCoefsExternal.iloc[0]["r"]
  cExternal = subCoefsExternal.iloc[0]["c"]
  fmtpExternal = fMTP(rExternal, cExternal, kExternal)
  
  # Subset participant's data by condition
  dataSubAction = dataAction[dataAction['ID']==sub]
  dataSubExternal = dataExternal[dataExternal['ID']==sub]
  
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
    mean_state_action['RT'] = temp2RT(mean_state_action.iloc[0]["prep"], subCoefsAction.iloc[0]["a SSE Action"], subCoefsAction.iloc[0]["b SSE Action"])
    
    # Predictions for external
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
          
    state_discr_external=state_discr_external[1:] # first trial has no preparation
    mean_state_external=state_discr_external.groupby(['FP'])[["prep"]].mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.iloc[0]["prep"], subCoefsExternal.iloc[0]["a SSE External"], subCoefsExternal.iloc[0]["b SSE External"])
          
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
  
  
  # Rename foreperiod columns in data DFs
  dataSubAction = dataSubAction.rename(columns = {"foreperiod":"FP"})
  dataSubExternal = dataSubExternal.rename(columns = {"foreperiod":"FP"})
  
  # Plot
  ax[iSub].plot(simAction.FP, simAction.RT, ".--", color = "blue")
  ax[iSub].scatter(dataSubAction.FP, dataSubAction.RT, marker = ".", color = "blue")
  ax[iSub].plot(simExternal.FP, simExternal.RT, ".--", color = "orange")
  ax[iSub].scatter(dataSubExternal.FP, dataSubExternal.RT, marker = ".", color = "orange")
  ax[iSub].set_title(sub)


plt.show()

figname = "./Modeling/fMTP/Fitting/Plots/con_FP_sse_fits_r.png"
plt.savefig(figname, format = "png", bbox_inches = "tight")

#======================= Compare parameters between conditions ======================
# Build dataset in long format for plotting
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop = True)
coefsAction = coefsAction[["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action"]]
coefsAction["condition"] = pd.Series(["action"] * len(coefsAction))
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop = True)
coefsExternal = coefsExternal[["ID", "k", "r", "c", "a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"]]
coefsExternal["condition"] = pd.Series(["external"] * len(coefsExternal))

coefsLong = pd.DataFrame(pd.concat([coefsAction.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                   "a SSE Action":"aSSE", "b SSE Action":"bSSE", 
                                                                   "SSEAction":"SSE",
                                                                   "R2_SSE_Action":"R2_SSE"}), 
                                   coefsExternal.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                     "a SSE External":"aSSE", 
                                                                     "b SSE External":"bSSE", 
                                                                     "SSEExternal": "SSE",
                                                                     "R2_SSE_External":"R2_SSE"})], 
                                  axis = 0), columns = ["ID", "condition", "k", "r", "c", "aSSE", "bSSE", "SSE", "R2_SSE"]).reset_index(drop = True)


# Scatter plot
paramFig, ax = plt.subplots(1,1)

sns.stripplot(x = "condition", y = "r", data = coefsLong, jitter = 0.1, orient = "v")

plt.show()

# Histogram
paramFig, ax = plt.subplots(1,1)

sns.histplot(x = "r", hue = "condition", data = coefsLong)

plt.show()

# Paired samples t-test
stats.ttest_rel(coefsAction.r, coefsExternal.r)


#===============================================================================#
#============================ 4. k and r fixed ==================================
#===============================================================================#

# Read empirical data and dataframe with fit results
subFitResults = pd.read_csv("./Modeling/fMTP/Fitting/subFitResults_c.csv")
empData = pd.read_csv("./Analysis/summaryData2.csv")
empData = empData.rename(columns = {"meanRT":"RT"})
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

# Separate datasets by condition 
dataAction = empData.loc[empData.condition == "action"].reset_index(drop=True)
dataExternal = empData.loc[empData.condition == "external"].reset_index(drop=True)

# Constants for simulations in both conditions
FPs = np.array([1.0, 2.8])
startVals=[4., 375.]

# Find minimum SSE for each condition
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop=True)
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop=True)

# Prepare panels for plotting
fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(4,6,sharey=False)
ax = ax.ravel()

titles = coefsAction.ID.tolist()

for iSub, sub in enumerate(titles):
  
  # Parameters for each condition for subject
  subCoefsAction = coefsAction[coefsAction.ID == sub]
  subCoefsExternal = coefsExternal[coefsExternal.ID == sub]
  
  kAction = subCoefsAction.iloc[0]["k"]
  rAction = subCoefsAction.iloc[0]["r"]
  cAction = subCoefsAction.iloc[0]["c"]
  fmtpAction = fMTP(rAction, cAction, kAction)
  
  kExternal = subCoefsExternal.iloc[0]["k"]
  rExternal = subCoefsExternal.iloc[0]["r"]
  cExternal = subCoefsExternal.iloc[0]["c"]
  fmtpExternal = fMTP(rExternal, cExternal, kExternal)
  
  # Subset participant's data by condition
  dataSubAction = dataAction[dataAction['ID']==sub]
  dataSubExternal = dataExternal[dataExternal['ID']==sub]
  
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
    mean_state_action['RT'] = temp2RT(mean_state_action.iloc[0]["prep"], subCoefsAction.iloc[0]["a SSE Action"], subCoefsAction.iloc[0]["b SSE Action"])
    
    # Predictions for external
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
          
    state_discr_external=state_discr_external[1:] # first trial has no preparation
    mean_state_external=state_discr_external.groupby(['FP'])[["prep"]].mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.iloc[0]["prep"], subCoefsExternal.iloc[0]["a SSE External"], subCoefsExternal.iloc[0]["b SSE External"])
          
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
  
  
  # Rename foreperiod columns in data DFs
  dataSubAction = dataSubAction.rename(columns = {"foreperiod":"FP"})
  dataSubExternal = dataSubExternal.rename(columns = {"foreperiod":"FP"})
  
  # Plot
  ax[iSub].plot(simAction.FP, simAction.RT, ".--", color = "blue")
  ax[iSub].scatter(dataSubAction.FP, dataSubAction.RT, marker = ".", color = "blue")
  ax[iSub].plot(simExternal.FP, simExternal.RT, ".--", color = "orange")
  ax[iSub].scatter(dataSubExternal.FP, dataSubExternal.RT, marker = ".", color = "orange")
  ax[iSub].set_title(sub)


plt.show()

figname = "./Modeling/fMTP/Fitting/Plots/con_FP_sse_fits_c.png"
plt.savefig(figname, format = "png", bbox_inches = "tight")

#======================= Compare parameters between conditions ======================
# Build dataset in long format for plotting
coefsAction = subFitResults.loc[subFitResults.groupby("ID").SSEAction.idxmin()].reset_index(drop = True)
coefsAction = coefsAction[["ID", "k", "r", "c", "a SSE Action", "b SSE Action", "SSEAction", "R2_SSE_Action"]]
coefsAction["condition"] = pd.Series(["action"] * len(coefsAction))
coefsExternal = subFitResults.loc[subFitResults.groupby("ID").SSEExternal.idxmin()].reset_index(drop = True)
coefsExternal = coefsExternal[["ID", "k", "r", "c", "a SSE External", "b SSE External", "SSEExternal", "R2_SSE_External"]]
coefsExternal["condition"] = pd.Series(["external"] * len(coefsExternal))

coefsLong = pd.DataFrame(pd.concat([coefsAction.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                   "a SSE Action":"aSSE", "b SSE Action":"bSSE", 
                                                                   "SSEAction":"SSE",
                                                                   "R2_SSE_Action":"R2_SSE"}), 
                                   coefsExternal.rename(columns = {"ID":"ID", "condition":"condition", "k":"k", "r":"r", "c":"c", 
                                                                     "a SSE External":"aSSE", 
                                                                     "b SSE External":"bSSE", 
                                                                     "SSEExternal": "SSE",
                                                                     "R2_SSE_External":"R2_SSE"})], 
                                  axis = 0), columns = ["ID", "condition", "k", "r", "c", "aSSE", "bSSE", "SSE", "R2_SSE"]).reset_index(drop = True)


# Scatter plot
paramFig, ax = plt.subplots(1,1)

sns.stripplot(x = "condition", y = "c", data = coefsLong, jitter = 0.1, orient = "v")

plt.show()

# Histogram
paramFig, ax = plt.subplots(1,1)

sns.histplot(x = "c", hue = "condition", data = coefsLong)

plt.show()

# Paired samples t-test
stats.ttest_rel(coefsAction.c, coefsExternal.c)

