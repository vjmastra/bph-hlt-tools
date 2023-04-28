from CRABClient.UserUtilities import config
config = config()

name = 'MuMuTrigger'
num = '_ZB'
era = '_362616'
tag = '_19feb23'

config.General.requestName = name + num + era + tag

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = '../test/PsikaonRootupler.py'

config.Data.inputDataset = '/ZeroBias/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/BcToJpsiPi_TuneCP5_14TeV_pythia8_evtgen/Run3Summer21MiniAOD-Pilot_120X_mcRun3_2021_realistic_v5_ext1-v2/MINIAODSIM'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = name + num + era + tag

# These values only make sense for processing data
#    Select input data based on a lumi mask
config.Data.lumiMask = './jsonFolder/Cert_Collisions2022_355100_362760_Golden.json'

# Select input data based on run-ranges
config.Data.runRange = '362616'

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_IT_Bari'
