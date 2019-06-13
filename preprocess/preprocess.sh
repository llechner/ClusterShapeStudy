# customized 10023.0
mkdir logs
cmsDriver.py MinBias_13TeV_pythia8_TuneCUETP8M1_cfi --conditions auto:phase1_2017_realistic -n 100 --era Run2_2017 --eventcontent FEVTDEBUG --relval 90000,100 -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:@relval2017 --datatier GEN-SIM-RAW --beamspot Realistic25ns13TeVEarly2017Collision --geometry DB:Extended --fileout file:step1.root > logs/step1.log  2>&1
cmsDriver.py step2  --conditions auto:phase1_2017_realistic -s RAW2DIGI,RECO:reconstruction_trackingOnly --datatier RECO -n 100 --geometry DB:Extended --era Run2_2017 --eventcontent FEVTDEBUG --filein  file:step1.root  --fileout file:step2.root > logs/step2.log  2>&1
rm MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.py
rm step1.root
cmsRun cfg/clusterShape_cfg.py > logs/clusterShape.log 2>&1
rm step2_RAW2DIGI_RECO.py
rm step2.root
