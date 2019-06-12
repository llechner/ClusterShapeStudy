import os
from CRABClient.UserUtilities import getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()

# Set variables
production_label    = "ClusterShapeStudy"
dataset             = "/Simu_Run3GT_test2/gbourgat-TTbar_step3_Run3GT_test3_RECOSIMoutput-09563c05d516fcc0d0acd0f72375719c/USER"
cfgFile             = "cfg/clusterShape_cfg.py"
unitsPerJob         = 10
publish             = False
runOnNonValid       = False
#totalUnits = 

config.section_("General")
config.General.transferLogs = True
config.General.requestName = production_label
config.General.workArea = 'crab_' + production_label

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = cfgFile
config.JobType.outputFiles = [production_label+'.root']

config.section_("Data")
config.Data.inputDataset = dataset

#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = False
config.Data.allowNonValidInputDataset = runOnNonValid

config.Data.outLFNDirBase = '/store/user/%s/tracking/%s/' % (getUsernameFromSiteDB(), production_label)
config.Data.publication = publish
config.Data.unitsPerJob = unitsPerJob
#if "CRAB_TOTAL_UNITS" in os.environ: config.Data.totalUnits = int(totalUnits)#8

config.section_("Site")
#config.Site.blacklist = ['T2_US_Purdue', 'T2_US_Nebraska', 'T2_US_MIT', 'T2_US_Caltech']
#config.Site.whitelist = ["T2_FR_IPHC"]
config.Site.storageSite = 'T2_AT_Vienna'

config.section_("User")

config.section_("Debug")




