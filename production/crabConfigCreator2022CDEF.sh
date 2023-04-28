#!/bin/bash

##
## just use : ./crabConfigCreator2022.sh
## works for 2022 C Dv1 Dv2 E F 
## change tag in the code
## add to pset two config files: one for 2022C-F
##

echo "************************ BEGIN ******************************"

PyFile=crabConfig_Data.py

#k = index of analyzer in psets
k=1

#i = C, Dv1, Dv2, E, F
for (( i=0; i<5; i++ ))
do

#j = 0, 1, 2, 3, 4, 5, 6, 7
for (( j=0; j<8; j++ ))
do

rm ${PyFile} 
echo "${PyFile} was deleted"
echo " "

echo "i = ${i}"
echo "j = ${j}"
echo " "

cat>> ${PyFile} <<pyFile

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

num = str(${j})

datasetnames = [
'/ParkingDoubleMuonLowMass' + num + '/Run2022C-PromptReco-v1/MINIAOD',
'/ParkingDoubleMuonLowMass' + num + '/Run2022D-PromptReco-v1/MINIAOD',
'/ParkingDoubleMuonLowMass' + num + '/Run2022D-PromptReco-v2/MINIAOD',
'/ParkingDoubleMuonLowMass' + num + '/Run2022E-PromptReco-v1/MINIAOD',
'/ParkingDoubleMuonLowMass' + num + '/Run2022F-PromptReco-v1/MINIAOD',
]

eras = [
'2022C',
'2022Dv1',
'2022Dv2',
'2022E',
'2022F',
]

pset = psets[$k]
tagname = tagnames[$k]
tag = '16feb23'

era = eras[$i]
dataset = datasetnames[$i]

config.General.requestName = tagname + '_' + era + '_DMLM' + num + '_' + tag

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = pset

config.Data.inputDataset = dataset
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 180
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

pyFile

echo "${PyFile} was created for ${i} and ${j}"
echo " "

crab submit ${PyFile} 

echo " "
echo "Crab task was submitted for ${i} and ${j}"
echo " "

done

done

echo "*************** END *****************"
