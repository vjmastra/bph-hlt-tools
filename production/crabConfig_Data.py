
from CRABClient.UserUtilities import config
config = config()

psets = [
'../test/PsikaonRootupler.py',
'../test/JPsiPiRootupler.py',
]

tagnames = [
'MuMuTrigger',
'Bc',
]

jsonfile = './jsonFolder/Cert_Collisions2022_355100_362760_Golden.json'

num = str(7)

datasetnames = [
'/ParkingDoubleMuonLowMass' + num + '/Run2022G-PromptReco-v1/MINIAOD',
]

eras = [
'2022G',
]

pset = psets[0]
tagname = tagnames[0]
tag = '19feb23'

era = eras[0]
dataset = datasetnames[0]

config.General.requestName = tagname + '_' + era + '_DMLM' + num + '_' + tag

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = pset

config.Data.inputDataset = dataset
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = tagname + '_' + era + '_DMLM' + num + '_' + tag

# These values only make sense for processing data
#    Select input data based on a lumi mask
config.Data.lumiMask = jsonfile

# Select input data based on run-ranges
config.Data.runRange = ''

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_IT_Bari'

