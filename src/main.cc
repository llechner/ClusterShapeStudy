#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <experimental/filesystem>
#include <dirent.h>
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <sys/types.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TRandom3.h"

#include "ClusterShape/interface/Estimator.h"
#include "ClusterShape/interface/ClusterShape.h"
//#include "DataFormats/Common/interface/DetSet.h"
//#include "DataFormats/Common/interface/DetSetVector.h"
//#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
//#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
//#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "Geometry/Records/interface/TrackerTopologyRcd.h"


bool PrintClusters = false;
bool PrintThresholds = false;
const double ClusterThreshold=40; // based on observation, make that flexible?
const int MaxClusterThreshold=70; // based on observation, make that flexible?
//const double ClusterThreshold=60.; // based on observation, make that flexible?


double x[14][3]={
    {0.848,0.063,0.013}, //TIB //IB1
    {0.878,0.053,0.008}, //IB2
    {0.808,0.079,0.017}, //TOB //OB2
    {0.769,0.093,0.023}, // OB1
    {0.857,0.0624,0.00911}, //TID //W1a
    {0.886,0.0502,0.00676}, //W2a
    {0.898,0.0496,0.00117}, //W3a
    {0.883,0.0528,0.00585}, //TEC //W1b
    {0.894,0.0489,0.0039}, //W2b
    {0.861,0.059,0.0104}, //W3b
    {0.888,0.0546,0.0013}, //W4
    {0.8,0.0804,0.0198}, //W5
    {0.807,0.0797,0.0169}, //W6
    {0.788,0.0913,0.0146}, //W7
};

double s[14][3]={
    {0.017,0.008,0.003}, //TIB //IB1
    {0.021,0.006,0.003}, //IB2
    {0.020,0.006,0.004}, //TOB //OB2
    {0.021,0.006,0.004}, //OB1
    {0.021,0.006,0.003}, //TID // W1a//s0=0.025 * x0, s1=0.1 * x1, s2=0.35 * x2 
    {0.022,0.005,0.002}, //W2a
    {0.022,0.005,0.0005}, //W3a
    {0.022,0.005,0.002}, //TEC //W1b
    {0.022,0.005,0.001}, //W2b
    {0.021,0.006,0.003}, //W3b
    {0.022,0.005,0.0005}, //W4
    {0.02,0.008,0.006}, //W5
    {0.02,0.008,0.005}, //W6
    {0.019,0.008,0.005}, //W7
};

const int TIBid=3;
const int TIDid=4;
const int TOBid=5;
const int TECid=6;

using namespace std;

int main(){

    //------------------------------------ //Ouverture fichier+SetBranch

    // const char *path="/eos/cms/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR17_Aag/";
    //string filename="/eos/cms/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR17_Aag/calibTree_304777_25.root";
    const string filename="data/input_v5.root";
    const string outfilename="data/results.root";

    const string plotDir="/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_datav4_v4_simhits2p_TEC/";
//    const string plotDir="/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_datav4_v4_simhits2p_TEC_dz005To04/";
//    const string plotDir="/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_datav4_v4_simhits2p_mesons_TOB/";
//    const string plotDir="/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_datav4_v4_simhits2p_TEC_dzGT04/";
    const string subPlotDir=plotDir+"/clusters/";

    try{
        std::experimental::filesystem::create_directory(plotDir);
    } catch (...){
    }

    try{
        std::experimental::filesystem::create_directory(subPlotDir);
    } catch (...){
    }

    try{
        std::experimental::filesystem::copy("data/index.php",plotDir.c_str());
    } catch (...){
    }

    try{
        std::experimental::filesystem::copy("data/index.php",subPlotDir.c_str());
    } catch (...){
    }

    TFile* ifile=TFile::Open(filename.c_str());
    TTree* tree=(TTree*) ifile->Get("testTree/tree");

    ostringstream title;

    std::vector<unsigned int>*GCStripIdx   = 0;
    std::vector<double>* GCAmplitudes      = 0;
    std::vector<unsigned int>*GCclusterIdx = 0;
    std::vector<double>* GCCharge          = 0;
    std::vector<unsigned int>*GCWidth      = 0;
    std::vector<double>* GCLocalTrackPhi   = 0;
    std::vector<double>* GCLocalTrackTheta = 0;
    std::vector<int>*GCSubdetid            = 0;
    std::vector<int>*GCLayerwheel          = 0;
    std::vector<unsigned int>*GCDetid      = 0;
    std::vector<double>* GCLocalpitch      = 0;
    std::vector<double>* GCSensorThickness = 0;
    std::vector<bool>*GCSaturation         = 0;
    std::vector<bool>*GCOverlapping        = 0;
    std::vector<bool>*GCFarfromedge        = 0;
    std::vector<double>* GCPath            = 0;

    std::vector<unsigned int>*tsosstrip = 0;
    std::vector<unsigned int>*tsostrackindex = 0;
    std::vector<unsigned int>*tsostrackmulti = 0;
    std::vector<unsigned int>*tsosclusterIdx = 0;
    std::vector<unsigned int>*tsosExpCW = 0;

    std::vector<unsigned int>*clusterStripIdx   = 0;
    std::vector<double>* clusterAmplitudes      = 0;
    std::vector<unsigned int>*clusterclusterIdx = 0;
    std::vector<double>* clusterCharge          = 0;
    std::vector<double>* clusterclusterCharge          = 0;
    std::vector<unsigned int>*clusterWidth      = 0;
    std::vector<int>*clusterSubdetid            = 0;
    std::vector<int>*clusterLayerwheel          = 0;
    std::vector<unsigned int>*clusterDetid      = 0;
    std::vector<double>* clusterLocalpitch      = 0;
    std::vector<double>* clusterSensorThickness = 0;
    std::vector<bool>*clusterSaturation         = 0;
    std::vector<bool>*clusterOverlapping        = 0;


    std::vector<unsigned int>*simstrip = 0;
    std::vector<double>*simtime1             = 0;
    std::vector<double>*simtime2             = 0;
    std::vector<double>*simtime3             = 0;
    std::vector<double>*simtime4             = 0;
    std::vector<unsigned int>*simhits       = 0;
    std::vector<double>* simLocalTrackPhi   = 0;
    std::vector<double>* simLocalTrackTheta = 0;

    std::vector<int>*simparticle1      = 0;
    std::vector<int>*simparticle2      = 0;
    std::vector<int>*simparticle3      = 0;
    std::vector<int>*simparticle4      = 0;

    std::vector<double>*simenergyloss1      = 0;
    std::vector<double>*simenergyloss2      = 0;
    std::vector<double>*simenergyloss3      = 0;
    std::vector<double>*simenergyloss4      = 0;

    std::vector<double>* simLocalx1   = 0;
    std::vector<double>* simLocaly1   = 0;
    std::vector<double>* simLocalz1   = 0;
    std::vector<double>* simLocalx2   = 0;
    std::vector<double>* simLocaly2   = 0;
    std::vector<double>* simLocalz2   = 0;

    std::vector<double>* simLocalx3   = 0;
    std::vector<double>* simLocaly3   = 0;
    std::vector<double>* simLocalz3   = 0;
    std::vector<double>* simLocalx4   = 0;
    std::vector<double>* simLocaly4   = 0;
    std::vector<double>* simLocalz4   = 0;

    tree->SetBranchAddress("GainCalibrationstripidx",&GCStripIdx);
    tree->SetBranchAddress("GainCalibrationamplitude",&GCAmplitudes);
    tree->SetBranchAddress("GainCalibrationclusteridx",&GCclusterIdx);
    tree->SetBranchAddress("GainCalibrationcharge",&GCCharge);
    tree->SetBranchAddress("GainCalibrationnstrips",&GCWidth);
    tree->SetBranchAddress("GainCalibrationlocalphi",&GCLocalTrackPhi);
    tree->SetBranchAddress("GainCalibrationlocaltheta",&GCLocalTrackTheta);
    tree->SetBranchAddress("GainCalibrationsubdetid",&GCSubdetid);
    tree->SetBranchAddress("GainCalibrationlayerwheel",&GCLayerwheel);
    tree->SetBranchAddress("GainCalibrationrawid",&GCDetid);
    tree->SetBranchAddress("GainCalibrationlocalpitch",&GCLocalpitch);
    tree->SetBranchAddress("GainCalibrationthickness",&GCSensorThickness);
    tree->SetBranchAddress("GainCalibrationsaturation",&GCSaturation);
    tree->SetBranchAddress("GainCalibrationoverlapping",&GCOverlapping);
    tree->SetBranchAddress("GainCalibrationfarfromedge",&GCFarfromedge);
    tree->SetBranchAddress("GainCalibrationpath",&GCPath);

    tree->SetBranchAddress("tsosstrip",&tsosstrip);
    tree->SetBranchAddress("tsostrackindex",&tsostrackindex);
    tree->SetBranchAddress("tsostrackmulti",&tsostrackmulti);
    tree->SetBranchAddress("tsosclusterIdx",&tsosclusterIdx);
    tree->SetBranchAddress("tsoscovered",&tsosExpCW);

    tree->SetBranchAddress("clusterstripidx",&clusterStripIdx);
    tree->SetBranchAddress("clusteramplitude",&clusterAmplitudes);
    tree->SetBranchAddress("clusterclusteridx",&clusterclusterIdx);
    tree->SetBranchAddress("clusterclustercharge",&clusterclusterCharge);
    tree->SetBranchAddress("clustercharge",&clusterCharge);
    tree->SetBranchAddress("clusterwidth",&clusterWidth);
    tree->SetBranchAddress("clustersubdetid",&clusterSubdetid);
    tree->SetBranchAddress("clusterlayerwheel",&clusterLayerwheel);
    tree->SetBranchAddress("clusterdetid",&clusterDetid);
    tree->SetBranchAddress("clusterlocalpitch",&clusterLocalpitch);
    tree->SetBranchAddress("clusterthickness",&clusterSensorThickness);
    tree->SetBranchAddress("clustersaturation",&clusterSaturation);
    tree->SetBranchAddress("clusteroverlapping",&clusterOverlapping);

    tree->SetBranchAddress("simstrip",&simstrip);
    tree->SetBranchAddress("simtime",&simtime1);
    tree->SetBranchAddress("simtime2",&simtime2);
    tree->SetBranchAddress("simtime3",&simtime3);
    tree->SetBranchAddress("simtime4",&simtime4);
    tree->SetBranchAddress("simhits",&simhits);
    tree->SetBranchAddress("simlocalphi",&simLocalTrackPhi);
    tree->SetBranchAddress("simlocaltheta",&simLocalTrackTheta);

    tree->SetBranchAddress("simparticle",&simparticle1);
    tree->SetBranchAddress("simparticle2",&simparticle2);
    tree->SetBranchAddress("simparticle3",&simparticle3);
    tree->SetBranchAddress("simparticle4",&simparticle4);

    tree->SetBranchAddress("simenergyloss",&simenergyloss1);
    tree->SetBranchAddress("simenergyloss2",&simenergyloss2);
    tree->SetBranchAddress("simenergyloss3",&simenergyloss3);
    tree->SetBranchAddress("simenergyloss4",&simenergyloss4);

    tree->SetBranchAddress("simlocalx",&simLocalx1);
    tree->SetBranchAddress("simlocaly",&simLocaly1);
    tree->SetBranchAddress("simlocalz",&simLocalz1);
    tree->SetBranchAddress("simlocalx2",&simLocalx2);
    tree->SetBranchAddress("simlocaly2",&simLocaly2);
    tree->SetBranchAddress("simlocalz2",&simLocalz2);

    tree->SetBranchAddress("simlocalx3",&simLocalx3);
    tree->SetBranchAddress("simlocaly3",&simLocaly3);
    tree->SetBranchAddress("simlocalz3",&simLocalz3);
    tree->SetBranchAddress("simlocalx4",&simLocalx4);
    tree->SetBranchAddress("simlocaly4",&simLocaly4);
    tree->SetBranchAddress("simlocalz4",&simLocalz4);

    //------------------------------------ //Declarations
    TCanvas* c1=new TCanvas("c1","c1");
    TCanvas* c2=new TCanvas("c2","c2");
    TCanvas* c3=new TCanvas("c3","c3");
    TCanvas* c4=new TCanvas("c4","c4",1200,800);
    TCanvas* c5=new TCanvas("c5","c5");
    TCanvas* c6=new TCanvas("c6","c6");
    TCanvas* c7=new TCanvas("c7","c7");
    TCanvas* c8=new TCanvas("c8","c8");
    TCanvas* c9=new TCanvas("c9","c9");
    TCanvas* c10=new TCanvas("c10","c10");
    TCanvas* c11=new TCanvas("c11","c11");
    TCanvas* c12=new TCanvas("c12","c12");
    TCanvas* c13=new TCanvas("c13","c13");
    TCanvas* c14=new TCanvas("c14","c14");
    TCanvas* c15=new TCanvas("c15","c15");
    TCanvas* c16=new TCanvas("c16","c16");
    TCanvas* c17=new TCanvas("c17","c17");
    TCanvas* c18=new TCanvas("c18","c18");
    TCanvas* c19=new TCanvas("c19","c19");
    TCanvas* c20=new TCanvas("c20","c20");
    TCanvas* c21=new TCanvas("c21","c21");
    TCanvas* c22=new TCanvas("c22","c22");
    TCanvas* c23=new TCanvas("c23","c23");
    TCanvas* c24=new TCanvas("c24","c24");
    TCanvas* c25=new TCanvas("c25","c25");
    TCanvas* c26=new TCanvas("c26","c26");

    TCanvas* c27=new TCanvas("c27","c27");
    TCanvas* c28=new TCanvas("c28","c28");
    TCanvas* c29=new TCanvas("c29","c29");
    TCanvas* c30=new TCanvas("c30","c30");
    TCanvas* c31=new TCanvas("c31","c31");
    TCanvas* c32=new TCanvas("c32","c32");
    TCanvas* c33=new TCanvas("c33","c33");
    TCanvas* c34=new TCanvas("c34","c34");
    TCanvas* c35=new TCanvas("c35","c35");
    TCanvas* c36=new TCanvas("c36","c36");
    TCanvas* c37=new TCanvas("c37","c37");
    TCanvas* c38=new TCanvas("c38","c38");
    TCanvas* c39=new TCanvas("c39","c39");
    TLegend* leg1=new TLegend(0.65,0.65,0.9,0.9);
    TLegend* leg2=new TLegend(0.65,0.65,0.9,0.9);
    TLegend* leg3=new TLegend(0.65,0.65,0.9,0.9);
    TLegend* leg4=new TLegend(0.65,0.65,0.9,0.9);
    TLegend* leg5=new TLegend(0.65,0.65,0.9,0.9);
    TLegend* leg6=new TLegend(0.55,0.75,0.9,0.9);
    TLegend* leg7=new TLegend(0.65,0.75,0.9,0.9);
    TLegend* leg8=new TLegend(0.65,0.75,0.9,0.9);
    TLegend* leg9=new TLegend(0.65,0.75,0.9,0.9);
    TLegend* leg10=new TLegend(0.65,0.65,0.9,0.9);

    TRandom3* r=new TRandom3();

    TH1F* hdEdx=new TH1F("hdEdx","",200,0,1000);

    double range[]={0,1.001,2.001,3.001,4.001,5.001,6.001,7.001,8.001};

    TH1F* hnTracks=new TH1F("hnTracks","",8,0,8);

    TH1F* hsubdetratio=new TH1F("hsubdetratio","",4,3,7);
    TH1F* hsubdet=new TH1F("hsubdet","",4,3,7);
    TH1F* hsubdet1=new TH1F("hsubdet1","",4,3,7);
    TH1F* hsubdet2p=new TH1F("hsubdet2p","",4,3,7);
    TH1F* hsubdetNorm1=new TH1F("hsubdetNorm1","",4,3,7);
    TH1F* hsubdetNorm2p=new TH1F("hsubdetNorm2p","",4,3,7);

    TH1F* hdistOverThick2=new TH1F("hdistOverThick2","",50,0,2);
    TH1F* hdistOverThick3=new TH1F("hdistOverThick3","",50,0,2);
    TH1F* hdistOverThick4=new TH1F("hdistOverThick4","",50,0,2);

    TH1F* hzdistOverThick2=new TH1F("hzdistOverThick2","",50,0,1);
    TH1F* hzdistOverThick3=new TH1F("hzdistOverThick3","",50,0,1);
    TH1F* hzdistOverThick4=new TH1F("hzdistOverThick4","",50,0,1);

    TH1F* h2DdistOverThick2=new TH1F("h2DdistOverThick2","",50,0,1);
    TH1F* h2DdistOverThick3=new TH1F("h2DdistOverThick3","",50,0,1);
    TH1F* h2DdistOverThick4=new TH1F("h2DdistOverThick4","",50,0,1);

    TH1F* hdeltaT=new TH1F("hdeltaT","",50,0,.003);
    TH1F* hdeltaT13=new TH1F("hdeltaT13","",50,0,.003);
    TH1F* hdeltaT14=new TH1F("hdeltaT14","",50,0,.003);

    TH1F* hdeltaTdet=new TH1F("hdeltaTdet","",50,0,.0001);
    TH1F* hdeltaT13det=new TH1F("hdeltaT13det","",50,0,.0001);
    TH1F* hdeltaT14det=new TH1F("hdeltaT14det","",50,0,.0001);

    TH1F* hdeltaV=new TH1F("hdeltaV","",50,0,50);

    TH1F* h2Dlocaldist2=new TH1F("h2Dlocaldist2","",50,0,0.03);
    TH1F* h2Dlocaldist3=new TH1F("h2Dlocaldist3","",50,0,0.03);
    TH1F* h2Dlocaldist4=new TH1F("h2Dlocaldist4","",50,0,0.03);

    TH1F* hlocaldist2=new TH1F("hlocaldist2","",50,0,0.03);
    TH1F* hlocaldist3=new TH1F("hlocaldist3","",50,0,0.03);
    TH1F* hlocaldist4=new TH1F("hlocaldist4","",50,0,0.03);

    TH1F* hzlocaldist2=new TH1F("hzlocaldist2","",50,0,0.03);
    TH1F* hzlocaldist3=new TH1F("hzlocaldist3","",50,0,0.03);
    TH1F* hzlocaldist4=new TH1F("hzlocaldist4","",50,0,0.03);
    TH1F* hsimhits=new TH1F("hsimhits","",3,1,4);

    TH1F* hwheelratio=new TH1F("hwheelratio","",9,1,10);
    TH1F* hwheel=new TH1F("hwheel","",9,1,10);
    TH1F* hwheel1=new TH1F("hwheel1","",9,1,10);
    TH1F* hwheel2p=new TH1F("hwheel2p","",9,1,10);
    TH1F* hwheelNorm1=new TH1F("hwheelNorm1","",9,1,10);
    TH1F* hwheelNorm2p=new TH1F("hwheelNorm2p","",9,1,10);

    TH1F* hthresUnfold=new TH1F("hthresUnfold","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hthresrec=new TH1F("hthresRec","",MaxClusterThreshold,0,MaxClusterThreshold);

    TH1F* hObsCW=new TH1F("hObsCW","",8,range);
    TH1F* hExpCW=new TH1F("hExpCW","",8,range);
    // TH1F* hExpCW=new TH1F("hExpCW","",200,0,8);
    TH1F* hUnfoldCW=new TH1F("hUnfoldCW","",8,range);
    TH1F* hRecCW=new TH1F("hRecCW","",8,range);
    TH1F* hRecSelCW=new TH1F("hRecSelCW","",8,range);
    TH1F* hRecSelCWUncert_P=new TH1F("hRecSelCWUncert_P","",8,range);
    TH1F* hRecSelCWUncert_M=new TH1F("hRecSelCWUncert_M","",8,range);
    TH1F* hWeightedWithUncertaintiesCrossTalkCW=new TH1F("hWeightedWithUncertaintiesCrossTalkCW","",8,range);
    TH1F* hRecSelProbCW=new TH1F("hRecSelProbCW","",8,range);

    TH1F* hSimCharge=new TH1F("hSimCharge","",50,0,600);
    TH1F* hSim1Charge=new TH1F("hSim1Charge","",50,1,600);
    TH1F* hSim2Charge=new TH1F("hSim2Charge","",50,1,600);
    TH1F* hSim3Charge=new TH1F("hSim3Charge","",50,1,600);
    TH1F* hSim4Charge=new TH1F("hSim4Charge","",50,1,600);
    TH1F* hObsCharge=new TH1F("hObsCharge","",50,0,600);
    TH1F* hUnfoldCharge=new TH1F("hUnfoldCharge","",50,0,600);
    TH1F* hRecCharge=new TH1F("hRecCharge","",50,0,600);
    TH1F* hRecSelCharge=new TH1F("hRecSelCharge","",50,0,600);

    TH1F* hRelDiffObs=new TH1F("hRelDiffObs","",100,-2,6);
    TH1F* hRelDiffUnfold=new TH1F("hRelDiffUnfold","",100,-2,6);
    TH1F* hRelDiffRec=new TH1F("hRelDiffRec","",100,-2,6);
    TH1F* hRelDiffRecSel=new TH1F("hRelDiffRecSel","",100,-2,6);
    TH1F* hRelDiffWeight=new TH1F("hRelDiffWeight","",100,-2,6);

    TH1F* hMeanRelDiffObs=new TH1F("hMeanRelDiffObs","",100,-2,6);
    TH1F* hMeanRelDiffUnfold=new TH1F("hMeanRelDiffUnfold","",100,-2,6);
    TH1F* hMeanRelDiffRec=new TH1F("hMeanRelDiffRec","",100,-2,6);
    TH1F* hMeanRelDiffRecSel=new TH1F("hMeanRelDiffRecSel","",100,-2,6);
    TH1F* hMeanRelDiffWeight=new TH1F("hMeanRelDiffWeight","",100,-2,6);

    TH1F* hRelDiffThreshObs=new TH1F("hRelDiffThreshObs","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hRelDiffThreshUnfold=new TH1F("hRelDiffThreshUnfold","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hRelDiffThreshRec=new TH1F("hRelDiffThreshRec","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hRelDiffThreshRecSel=new TH1F("hRelDiffThreshRecSel","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hRelDiffThreshWeight=new TH1F("hRelDiffThreshWeight","",MaxClusterThreshold,0,MaxClusterThreshold);

    TH1F* hWidthThreshObs=new TH1F("hWidthThreshObs","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hWidthThreshUnfold=new TH1F("hWidthThreshUnfold","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hWidthThreshRec=new TH1F("hWidthThreshRec","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hWidthThreshRecSel=new TH1F("hWidthThreshRecSel","",MaxClusterThreshold,0,MaxClusterThreshold);
    TH1F* hWidthThreshWeight=new TH1F("hWidthThreshWeight","",MaxClusterThreshold,0,MaxClusterThreshold);

    TH1F* hRelDiff=new TH1F("hRelDiff","",100,-2,6);
    TH1F* hRelDiff1=new TH1F("hRelDiff1","",100,-2,6);
    TH1F* hRelDiff12=new TH1F("hRelDiff12","",100,-2,6);
    TH1F* hRelDiff2=new TH1F("hRelDiff2","",100,-2,6);
    TH1F* hRelDiff23=new TH1F("hRelDiff23","",100,-2,6);
    TH1F* hRelDiff3=new TH1F("hRelDiff3","",100,-2,6);

    // TH1F* hRelDiffWeighted=new TH1F("hRelDiffWeighted","",100,-0.5,1.1);
    TH1F* hRelDiff1Weighted=new TH1F("hRelDiff1Weighted","",100,-0.5,1.1);
    TH1F* hRelDiff12Weighted=new TH1F("hRelDiff12Weighted","",100,-0.5,1.1);
    TH1F* hRelDiff2Weighted=new TH1F("hRelDiff2Weighted","",100,-0.5,1.1);
    TH1F* hRelDiff23Weighted=new TH1F("hRelDiff23Weighted","",100,-0.5,1.1);
    TH1F* hRelDiff3Weighted=new TH1F("hRelDiff3Weighted","",100,-0.5,1.1);

    TH1F* hRelDiffCrossTalk_P=new TH1F("hRelDiffCrossTalk_P","",100,-0.5,1.1);
    TH1F* hRelDiffCrossTalk_M=new TH1F("hRelDiffCrossTalk_M","",100,-0.5,1.1);
    TH1F* hRelDiffCrossTalkWeighted=new TH1F("hRelDiffCrossTalkWeighted","",100,-0.5,1.1);
    TH1F* hRelDiffCrossTalkProb=new TH1F("hRelDiffCrossTalkProb","",100,-0.5,1.1);

    TH2F* hExpVsRec=new TH2F("hExpVsRec","",500,0,5,500,0,5);
    TH2F* hExpVsRecSel=new TH2F("hExpVsRecSel","",500,0,5,500,0,5);
    TH2F* hExpVsObs=new TH2F("hExpVsObs","",500,0,5,500,0,5);
    TH2F* hExpVsUnfold=new TH2F("hExpVsUnfold","",500,0,5,500,0,5);
    TH2F* hExpVsWeight=new TH2F("hExpVsWeight","",500,0,5,500,0,5);

    TProfile* pObsCluster=new TProfile("pObsCluster","Observed",16,0,15,-50,700);
    TProfile* pUnfoldCluster=new TProfile("pUnfoldCluster","Unfolded",16,0,15,-50,700);

    TProfile* pUnderCluster=new TProfile("pSubCluster","Sub-Cluster",16,0,15,-50,700);
    TProfile* pClusterSel=new TProfile("pClusterSel","Cleaned Cluster",16,0,15,-50,700);
    TProfile* ligneThreshold=new TProfile("","",1,0,1);
    TProfile* ligneQStrip=new TProfile("","",1,0,1);

    //Estimator Estim;

    double x0=0, x1=0, x2=0;
    double s0=0, s1=0;//, s2;

    unsigned int Idx=0;
    unsigned int cIdx=0;
    unsigned int tcIdx=0;

    unsigned int  k=0.;

    int binmax=0.;
    double  binx=0.;

    int subdet=0.;
    int ntracks=0.;
    int trkidx=0.;
    int shits=0.;

    int pdg1=-999;
    int pdg2=-999;
    int pdg3=-999;
    int pdg4=-999;

    double eLoss1=0.;
    double eLoss2=0.;
    double eLoss3=0.;
    double eLoss4=0.;
    double eLossTot=0.;

    double ObsCW=0.;
    double ExpCW=0.;
    double UnfoldCW=0.;
    double Unfold13CW=0.;
    double RecCW=0.;
    double RecSelCW=0.;
    double RecSelCWUncert_P=0.;
    double RecSelCWUncert_M=0.;
    double WeightedWithUncertaintiesCrossTalkCW=0.;
    double RecSelProbCW=0.;

    double ObsCharge=0.;
    double UnfoldCharge=0.;
    double RecCharge=0.;
    double RecSelCharge=0.;

    double RelDiffObs=0.;
    double RelDiffUnfold=0.;
    double RelDiffRec=0.;
    double RelDiffRecSel=0.;
    double RelDiffWeight=0.;
    double RelDiff=0.;
    double RelDiff1=0.;
    double RelDiff12=0.;
    double RelDiff2=0.;
    double RelDiff23=0.;
    double RelDiff3=0.;

    double RelDiffCrossTalk_P=0.;
    double RelDiffCrossTalk_M=0.;
    double RelDiffCrossTalkWeighted=0.;
    double RelDiffCrossTalkProb=0.;

    double wheel=0.;
    double sstrip=0.;
    double strip=0.;
    double QStrip=0.;
    double qstripfromtreatment=0.;

    double thick=0.;

    double deltaT12=0.;
    double deltaT13=0.;
    double deltaT14=0.;
    double deltaV=0.;

    double localdist2=0.;
    double localdist3=0.;
    double localdist4=0.;

    double zlocaldist2=0.;
    double zlocaldist3=0.;
    double zlocaldist4=0.;

    double localdist2d2=0.;
    double localdist2d3=0.;
    double localdist2d4=0.;

    double slx1=0.;
    double sly1=0.;
    double slz1=0.;
    double slx2=0.;
    double sly2=0.;
    double slz2=0.;
    double slx3=0.;
    double sly3=0.;
    double slz3=0.;
    double slx4=0.;
    double sly4=0.;
    double slz4=0.;

    std::vector<int> clusterIds;
    clusterIds.reserve(10000);
    std::vector<int> alltracks;
    alltracks.reserve(10000);

    std::vector<double> q;
    std::vector<double> QUnfold;
    std::vector<double> QUnfold_P;
    std::vector<double> QUnfold_M;
    std::vector<double> QUnfoldProb;

    std::vector<std::vector<double>> Cluster;
    std::vector<std::vector<double>> Cluster13;
    std::vector<double> UnderCluster;

    //------------------------------------

    //------------------------------------ RelUncertaintyCrossTalkNoise
    // double UncertQStrip1=0.581025; //1%
    // double UncertQStrip5=0.378116; //5%
    double UncertQStrip10=0.268006; //10%
    // double UncertQStrip20=0.139197; //20%
    double UncertQS=0.;
    //------------------------------------

//    edm::ESHandle<TrackerTopology> tTopoHandle;
//    const TrackerTopology* const tTopo = tTopoHandle.product();

    for (unsigned int th=0;th<=MaxClusterThreshold;th++){
      if(!PrintThresholds and th!=ClusterThreshold) continue;
      for (unsigned int i=0;i<tree->GetEntries();i++){
        tree->GetEntry(i);
        if(i%1000==0) cout<<"Event "<<i<<endl;

        for(unsigned int c=0;c<GCWidth->size();c++){

            if(GCSubdetid->at(c)==TIBid){
                 if(GCLayerwheel->at(c)<=2){
                     x0=x[0][0];
                     x1=x[0][1];
                     x2=x[0][2];
                     s0=s[0][0];
                     s1=s[0][1];
                 }
                 else {
                     x0=x[1][0];
                     x1=x[1][1];
                     x2=x[1][2];
                     s0=s[1][0];
                     s1=s[1][1];
                 }
             }

            if(GCSubdetid->at(c)==TIDid){
                if(GCLayerwheel->at(c)==1){
                    x0=x[4][0];
                    x1=x[4][1];
                    x2=x[4][2];
                }
                else if(GCLayerwheel->at(c)==2){
                    x0=x[5][0];
                    x1=x[5][1];
                    x2=x[5][2];
                }
                else if(GCLayerwheel->at(c)==3){
                    x0=x[6][0];
                    x1=x[6][1];
                    x2=x[6][2];
                }
            }

            if(GCSubdetid->at(c)==TOBid){
                if(GCLayerwheel->at(c)<=4){
                    x0=x[2][0];
                    x1=x[2][1];
                    x2=x[2][2];
                    s0=s[2][0];
                    s1=s[2][1];
                }
                else {
                    x0=x[3][0];
                    x1=x[3][1];
                    x2=x[3][2];
                    s0=s[3][0];
                    s1=s[3][1];
                }
            }

            if(GCSubdetid->at(c)==TECid){
                if(GCLayerwheel->at(c)==1){
                    x0=x[7][0];
                    x1=x[7][1];
                    x2=x[7][2];
                }
                else if(GCLayerwheel->at(c)==2){
                    x0=x[8][0];
                    x1=x[8][1];
                    x2=x[8][2];
                }
                else if(GCLayerwheel->at(c)==3){
                    x0=x[9][0];
                    x1=x[9][1];
                    x2=x[9][2];
                }
                else if(GCLayerwheel->at(c)==4){
                    x0=x[10][0];
                    x1=x[10][1];
                    x2=x[10][2];
                }
                else if(GCLayerwheel->at(c)==5){
                    x0=x[11][0];
                    x1=x[11][1];
                    x2=x[11][2];
                }
                else if(GCLayerwheel->at(c)==6){
                    x0=x[12][0];
                    x1=x[12][1];
                    x2=x[12][2];
                }
                else if(GCLayerwheel->at(c)==7){
                    x0=x[13][0];
                    x1=x[13][1];
                    x2=x[13][2];
                }
            }

            if(GCOverlapping->at(c)==0 && GCFarfromedge->at(c)==1 && GCSaturation->at(c)==0){

                tcIdx=tsosclusterIdx->at(c);
                strip=tsosstrip->at(c);
                ntracks=tsostrackmulti->at(c);
                trkidx=tsostrackindex->at(c);

                sstrip=simstrip->at(tcIdx);
                shits=simhits->at(tcIdx);
                subdet=clusterSubdetid->at(tcIdx);
                wheel=clusterLayerwheel->at(tcIdx);

                ObsCW=clusterWidth->at(tcIdx);
                ObsCharge=clusterclusterCharge->at(tcIdx); //Sum(q);

                pdg1 = simparticle1->at(tcIdx);
                pdg2 = simparticle2->at(tcIdx);
                pdg3 = simparticle3->at(tcIdx);
                pdg4 = simparticle4->at(tcIdx);

                // Only mesons in merged clusters
//                if( abs(pdg1)<30 || abs(pdg2)<30 || abs(pdg3)<30 || abs(pdg4)<30 ) continue;

//                cout<<"p1 "<<pdg1<<" dz1 "<<slz1<<" p2 "<<pdg2<<" dz2 "<<slz2<<" p3 "<<pdg3<<" dz3 "<<slz3<<" p4 "<<pdg4<<" dz4 "<<slz4<<endl;
//                if(shits != 1) continue;
                if(shits < 2) continue;
//                if(shits < 1) continue;

//                if (GCSubdetid->at(c) != TOBid) continue;
//                if (GCSubdetid->at(c) != TIBid) continue;
//                if (GCSubdetid->at(c) != TIDid) continue;
                if (GCSubdetid->at(c) != TECid) continue;
                //eLossTot = 0;

                if(shits > 0) {
                    eLoss1 = simenergyloss1->at(tcIdx)*1000000;
                    eLossTot = eLoss1;
                }
                if(shits > 1) {
                    eLoss2 = simenergyloss2->at(tcIdx)*1000000;
//                    eLossTot += eLoss2;
                }
                if(shits > 2) {
                    eLoss3 = simenergyloss3->at(tcIdx)*1000000;
//                    eLossTot += eLoss3;
                }
                if(shits > 3) {
                    eLoss4 = simenergyloss4->at(tcIdx)*1000000;
//                    eLossTot += eLoss4;
                }

//                cout<<"ELoss1 "<<eLoss1<<endl;
//                cout<<"ELoss2 "<<eLoss2<<endl;
//                cout<<"ObsCharge "<<ObsCharge<<endl;
//                cout<<endl;

                thick = GCSensorThickness->at(c);

                if(shits > 0) {
                    slx1 = simLocalx1->at(tcIdx);
                    sly1 = simLocaly1->at(tcIdx);
                    slz1 = simLocalz1->at(tcIdx);
                }
                if(shits > 1) {
                    slx2 = simLocalx2->at(tcIdx);
                    sly2 = simLocaly2->at(tcIdx);
                    slz2 = simLocalz2->at(tcIdx);
                    localdist2d2 = sqrt( pow((slx1-slx2), 2) + pow((sly1-sly2), 2) );
                    zlocaldist2 = abs(slz1-slz2);
                    localdist2 = sqrt( pow((slx1-slx2), 2) + pow((sly1-sly2), 2) + pow((slz1-slz2), 2) );
                    deltaT12 = simtime2->at(tcIdx) - simtime1->at(tcIdx);
//                    if(deltaT12<0.00001) cout<<deltaT12<<" "<<localdist2<<endl;
                    deltaV = thick / deltaT12;
                }
                if(shits > 2) {
                    slx3 = simLocalx3->at(tcIdx);
                    sly3 = simLocaly3->at(tcIdx);
                    slz3 = simLocalz3->at(tcIdx);
                    localdist2d3 = sqrt( pow((slx1-slx3), 2) + pow((sly1-sly3), 2) );
                    zlocaldist3 = abs(slz1-slz3);
                    localdist3 = sqrt( pow((slx1-slx3), 2) + pow((sly1-sly3), 2) + pow((slz1-slz3), 2) );
                    deltaT13 = simtime3->at(tcIdx) - simtime1->at(tcIdx);
                }
                if(shits > 3) {
                    slx4 = simLocalx4->at(tcIdx);
                    sly4 = simLocaly4->at(tcIdx);
                    slz4 = simLocalz4->at(tcIdx);
                    localdist2d4 = sqrt( pow((slx1-slx4), 2) + pow((sly1-sly4), 2) );
                    zlocaldist4 = abs(slz1-slz4);
                    localdist4 = sqrt( pow((slx1-slx4), 2) + pow((sly1-sly4), 2) + pow((slz1-slz4), 2) );
                    deltaT14 = simtime4->at(tcIdx) - simtime1->at(tcIdx);
                }

//                if(zlocaldist2/thick > 0.4 || zlocaldist2/thick < 0.05) continue;
//                if(zlocaldist2/thick < 0.4) continue;
//                if(zlocaldist2/thick > 0.05) continue;

                if(th==ClusterThreshold){
                     hsimhits->Fill(shits);
                    if(shits == 1){
                        hsubdet1->Fill(subdet);
                        hsubdet->Fill(subdet);
                        hwheel1->Fill(wheel);
                        hwheel->Fill(wheel);
                    } else if(shits > 1){
                        hsubdet->Fill(subdet);
                        hsubdet2p->Fill(subdet);
                        hsubdetratio->Fill(subdet);
                        hwheel->Fill(wheel);
                        hwheel2p->Fill(wheel);
                        hwheelratio->Fill(wheel);
                    }
                }


                for(unsigned int z=0;z<clusterStripIdx->size();z++){
                    if(clusterStripIdx->at(z) == clusterclusterIdx->at(tcIdx)){
                        q.push_back(clusterAmplitudes->at(z));
                    }
                }


//                if(zlocaldist2/thick>0.05) continue;

                double theta_L = 0.02;
                ExpCW= abs( tan(theta_L) + cos(GCLocalTrackPhi->at(c)) * tan(GCLocalTrackTheta->at(c)) ) * thick / GCLocalpitch->at(c);


                double s1b=UncertaintyCrossTalk(x0,s0,x1,x2);
                QUnfold=Unfold(q,x0,x1,x2);
//                UnfoldCharge=Sum(QUnfold);
                QUnfold_M=Unfold(q,x0+s0,x1+s1b,x2);
                QUnfold_P=Unfold(q,x0-s0,x1-s1b-2*s0,x2);






                double x0prob=r->Gaus(x0,2*s0);
                double x1prob=r->Gaus(x1,2*s1);
                double x2prob=(double)(1-2*x1-x0)/2.;

                QUnfoldProb=Unfold(q,x0prob,x1prob,x2prob);
                Cluster=Cutting(Threshold(QUnfold,th));
                UnfoldCharge=Sum(Threshold(QUnfold,th));

                if(Cluster.size()==0) cerr<<"Warning"<<"\tEvent "<<i<<"\tCluster"<<c<<endl;
                else UnderCluster=SelectionAfterCutting(Cluster);

                for(unsigned int i=0;i<Cluster.size();i++){
                    UnfoldCW += Cluster[i].size();
                }

                k=Cluster.size();
                if (k>1 && k<10) UnfoldCW+=k-1;

                Cluster13=Cutting(Threshold(QUnfold,13));

                for(unsigned int i=0;i<Cluster13.size();i++){
                    Unfold13CW += Cluster13[i].size();
                }

                k=Cluster13.size();
                if (k>1 && k<10) Unfold13CW+=k-1;

                QStrip=funcQStrip(UnderCluster,GCPath,c);

                std::vector<std::vector<double>> TotalClusterWidthVec=TotalClusterWidth(UnderCluster,QStrip);
                RecCW=TotalClusterWidthVec[1][0];
                RecCharge=TotalClusterWidthVec[1][1];

                UncertQS=UncertQStrip10;

                std::vector<double> TotalClusterWidthAfterTreatVec=TotalClusterWidthAfterTreatment(QUnfold,UncertQS,GCPath,c,th);
                RecSelCW=TotalClusterWidthAfterTreatVec[0];
                RecSelCharge=TotalClusterWidthAfterTreatVec[1];
                RecSelCWUncert_P=TotalClusterWidthAfterTreatment(QUnfold_P,UncertQS,GCPath,c,th)[0];
                RecSelCWUncert_M=TotalClusterWidthAfterTreatment(QUnfold_M,UncertQS,GCPath,c,th)[0];

                RelDiffObs=(double)(ObsCW-ExpCW)/ExpCW;
                RelDiffRec=(double)(RecCW-ExpCW)/ExpCW;
                RelDiffRecSel=(double)(RecSelCW-ExpCW)/ExpCW;
                RelDiffUnfold=(double)(UnfoldCW-ExpCW)/ExpCW;
                RelDiffWeight=(double)(WeightedWithUncertaintiesCrossTalkCW-ExpCW)/ExpCW;
                RelDiffCrossTalk_P=(double)(RecSelCWUncert_P-RecSelCW)/RecSelCW;
                RelDiffCrossTalk_M=(double)(RecSelCWUncert_M-RecSelCW)/RecSelCW;
                WeightedWithUncertaintiesCrossTalkCW=0.5*RecSelCW+0.25*(RecSelCWUncert_M+RecSelCWUncert_P);
                RelDiffCrossTalkProb=(double)(RecSelProbCW-ExpCW)/ExpCW;
                RelDiffCrossTalkWeighted=(double)(WeightedWithUncertaintiesCrossTalkCW-RecSelCW)/RecSelCW;

                if(RecSelCW==1){
                    RelDiff1=(double)(RecSelCW-ExpCW)/ExpCW;
                    hRelDiff1->Fill(RelDiff1);
                    hRelDiff1Weighted->Fill(RelDiffCrossTalkWeighted);
                }
                if(RecSelCW>1 && RecSelCW<2){
                    RelDiff12=(double)(RecSelCW-ExpCW)/ExpCW;
                    hRelDiff12->Fill(RelDiff12);
                    hRelDiff12Weighted->Fill(RelDiffCrossTalkWeighted);
                }
                if(RecSelCW==2){
                    RelDiff2=(double)(RecSelCW-ExpCW)/ExpCW;
                    hRelDiff2->Fill(RelDiff2);
                    hRelDiff2Weighted->Fill(RelDiffCrossTalkWeighted);
                }
                if(RecSelCW>2 && RecSelCW<3){
                    RelDiff23=(double)(RecSelCW-ExpCW)/ExpCW;
                    hRelDiff23->Fill(RelDiff23);
                    hRelDiff23Weighted->Fill(RelDiffCrossTalkWeighted);
                }
                if(RecSelCW>=3){
                    RelDiff3=(double)(RecSelCW-ExpCW)/ExpCW;
                    hRelDiff3->Fill(RelDiff3);
                    hRelDiff3Weighted->Fill(RelDiffCrossTalkWeighted);
                }

                vector<double> ClusterSel=ClusterizationWithNoise(UnderCluster,QStrip,UncertQS);

//                if(Unfold13CW != UnfoldCW && PrintClusters && Cluster.size()!=0){
//                    cout<<"13 "<<Unfold13CW<<" th "<<UnfoldCW<<" c "<<c<<endl;
                if(PrintClusters){

                    title.str("");
                    title.clear();
                    title << "Observed (SimHits = " << shits << ")";
                    pObsCluster->SetTitle( title.str().c_str() );
                    title.str("");
                    title.clear();
                    title << "Unfolded (SimHits = " << shits << ")";
                    pUnfoldCluster->SetTitle( title.str().c_str() );
                    title.str("");
                    title.clear();
                    title << "Sub-Cluster (SimHits = " << shits << ")";
                    pUnderCluster->SetTitle( title.str().c_str() );
                    title.str("");
                    title.clear();
                    title << "Cleaned Cluster (SimHits = " << shits << ")";
                    pClusterSel->SetTitle( title.str().c_str() );

                    pObsCluster->SetBins(q.size()+2,0,q.size()+2);
                    pUnfoldCluster->SetBins(QUnfold.size()+2,0,QUnfold.size()+2);
                    pUnderCluster->SetBins(UnderCluster.size()+2,0,UnderCluster.size()+2);
                    pClusterSel->SetBins(ClusterSel.size()+2,0,ClusterSel.size()+2);
                    ligneThreshold->SetBins(q.size()+2,0,q.size()+2);
                    ligneQStrip->SetBins(UnderCluster.size()+2,0,UnderCluster.size()+2);

                    for(unsigned int l=0;l<q.size();l++){
                        pObsCluster->Fill(l+1,q[l]);
                        pUnfoldCluster->Fill(l+1,QUnfold[l]);
                    }
                    for(unsigned int l=0;l<q.size()+2;l++) ligneThreshold->Fill(l,th);
                    pUnfoldCluster->Reset();
                    for(unsigned int l=0;l<QUnfold.size();l++){
                        if(QUnfold[l]<=th) QUnfold[l]=0.;
                        pUnfoldCluster->Fill(l+1,QUnfold[l]);
                    }
                    for(unsigned int l=0;l<UnderCluster.size();l++){
                        pUnderCluster->Fill(l+1,UnderCluster[l]);
                    }

                    for(unsigned int l=0;l<UnderCluster.size()+2;l++) ligneQStrip->Fill(l,qstripfromtreatment);
                    for(unsigned int l=0;l<ClusterSel.size();l++){
                        pClusterSel->Fill(l+1,ClusterSel[l]);
                    }


                    TLegend* legend=new TLegend(0.65,0.7,0.9,0.9);
                    c4->cd();
                    c4->Divide(2,2);
                    gStyle->SetOptStat(0);
                    c4->cd(1);
                    pObsCluster->GetXaxis()->SetTitle("Strip");
                    pObsCluster->GetYaxis()->SetTitle("Charge");
                    pObsCluster->GetYaxis()->SetTitleOffset(0.95);
                    pObsCluster->GetYaxis()->SetRangeUser(0,pUnfoldCluster->GetMaximum()+9);
                    pObsCluster->SetLineColor(4);
                    pObsCluster->SetLineWidth(2);

                    pObsCluster->GetXaxis()->SetLabelSize(0.045);
                    pObsCluster->GetYaxis()->SetLabelSize(0.045);
                    pObsCluster->GetXaxis()->SetTitleSize(0.045);
                    pObsCluster->GetYaxis()->SetTitleSize(0.045);

                    pObsCluster->Draw();
                    pUnfoldCluster->SetLineColor(2);
                    pUnfoldCluster->SetLineWidth(2);
                    pUnfoldCluster->Draw("SAME");
                    ligneThreshold->SetLineColor(1);
                    ligneThreshold->SetLineWidth(2);
                    ligneThreshold->SetLineStyle(7);
                    ligneThreshold->Draw("SAME");
                    legend->AddEntry(pObsCluster,"Observed","l");
                    legend->AddEntry(pUnfoldCluster,"Unfold","l");
                    legend->Draw("SAME");
                    c4->cd(2);
                    pUnfoldCluster->GetXaxis()->SetTitle("Strip");
                    pUnfoldCluster->GetYaxis()->SetTitle("Charge");
                    pUnfoldCluster->GetYaxis()->SetTitleOffset(0.95);

                    pUnfoldCluster->GetXaxis()->SetLabelSize(0.045);
                    pUnfoldCluster->GetYaxis()->SetLabelSize(0.045);
                    pUnfoldCluster->GetXaxis()->SetTitleSize(0.045);
                    pUnfoldCluster->GetYaxis()->SetTitleSize(0.045);

                    pUnfoldCluster->Draw();
                    c4->cd(3);
                    pUnderCluster->GetXaxis()->SetTitle("Strip");
                    pUnderCluster->GetYaxis()->SetTitle("Charge");
                    pUnderCluster->GetYaxis()->SetTitleOffset(0.95);

                    pUnderCluster->GetXaxis()->SetLabelSize(0.045);
                    pUnderCluster->GetYaxis()->SetLabelSize(0.045);
                    pUnderCluster->GetXaxis()->SetTitleSize(0.045);
                    pUnderCluster->GetYaxis()->SetTitleSize(0.045);

                    pUnderCluster->Draw();
                    ligneQStrip->SetLineColor(1);
                    ligneQStrip->SetLineWidth(2);
                    ligneQStrip->SetLineStyle(7);
                    ligneQStrip->Draw("SAME");
                    c4->cd(4);
                    pClusterSel->GetXaxis()->SetTitle("Strip");
                    pClusterSel->GetYaxis()->SetTitle("Charge");
                    pClusterSel->GetYaxis()->SetTitleOffset(0.95);

                    pClusterSel->GetXaxis()->SetLabelSize(0.045);
                    pClusterSel->GetYaxis()->SetLabelSize(0.045);
                    pClusterSel->GetXaxis()->SetTitleSize(0.045);
                    pClusterSel->GetYaxis()->SetTitleSize(0.045);

                    pClusterSel->Draw();
                    c4->SaveAs((subPlotDir+"chargeCluster_"+to_string(c)+".pdf").c_str());
                    c4->SaveAs((subPlotDir+"chargeCluster_"+to_string(c)+".png").c_str());
//                    getchar();
                    pObsCluster->Reset();
                    pUnfoldCluster->Reset();
                    pUnderCluster->Reset();
                    pClusterSel->Reset();
                    ligneThreshold->Reset();
                    ligneQStrip->Reset();
                    c4->Clear();

                }

                hMeanRelDiffObs->Fill(RelDiffObs);
                hMeanRelDiffUnfold->Fill(RelDiffUnfold);
                hMeanRelDiffRec->Fill(RelDiffRec);
                hMeanRelDiffRecSel->Fill(RelDiffRecSel);
                hMeanRelDiffWeight->Fill(RelDiffWeight);

                if(th==ClusterThreshold){

                    if(shits > 1){
                        hzlocaldist2->Fill(zlocaldist2);
                        h2Dlocaldist2->Fill(localdist2d2);
                        h2DdistOverThick2->Fill(localdist2d2/thick);
                        hlocaldist2->Fill(localdist2);
                        hdistOverThick2->Fill(localdist2/thick);
                        hzdistOverThick2->Fill(zlocaldist2/thick);
                        hdeltaT->Fill(abs(deltaT12));
                        hdeltaTdet->Fill(abs(deltaT12));
                        hdeltaV->Fill(deltaV);
                    }
                    if(shits > 2) {
                        hzlocaldist3->Fill(zlocaldist3);
                        h2Dlocaldist3->Fill(localdist2d3);
                        h2DdistOverThick3->Fill(localdist2d3/thick);
                        hlocaldist3->Fill(localdist3);
                        hdistOverThick3->Fill(localdist3/thick);
                        hzdistOverThick3->Fill(zlocaldist3/thick);
                        hdeltaT13->Fill(abs(deltaT13));
                        hdeltaT13det->Fill(abs(deltaT13));
                    }
                    if(shits > 3) {
                        hzlocaldist4->Fill(zlocaldist4);
                        h2Dlocaldist4->Fill(localdist2d4);
                        h2DdistOverThick4->Fill(localdist2d4/thick);
                        hlocaldist4->Fill(localdist4);
                        hdistOverThick4->Fill(localdist4/thick);
                        hzdistOverThick4->Fill(zlocaldist4/thick);
                        hdeltaT14->Fill(abs(deltaT14));
                        hdeltaT14det->Fill(abs(deltaT14));
                    }
                    //hsimhits->Fill(shits);

                    hnTracks->Fill(ntracks);
                    hObsCW->Fill(ObsCW);
                    hExpCW->Fill(ExpCW);
                    hUnfoldCW->Fill(UnfoldCW);
                    hRecCW->Fill(RecCW);
                    hRecSelCW->Fill(RecSelCW);
                    hRecSelCWUncert_P->Fill(RecSelCWUncert_P);
                    hRecSelCWUncert_M->Fill(RecSelCWUncert_M);
                    hRelDiff->Fill(RelDiff);
                    hRelDiffCrossTalk_P->Fill(RelDiffCrossTalk_P);
                    hRelDiffCrossTalk_M->Fill(RelDiffCrossTalk_M);
                    hRelDiffCrossTalkProb->Fill(RelDiffCrossTalkProb);
                    hExpVsRec->Fill(ExpCW,RecCW);
                    hExpVsObs->Fill(ExpCW,ObsCW);
                    hExpVsRecSel->Fill(ExpCW,RecSelCW);
                    hExpVsUnfold->Fill(ExpCW,UnfoldCW);
                    hExpVsWeight->Fill(ExpCW,WeightedWithUncertaintiesCrossTalkCW);
                    hWeightedWithUncertaintiesCrossTalkCW->Fill(WeightedWithUncertaintiesCrossTalkCW);
                    hRecSelProbCW->Fill(RecSelProbCW);
                    hObsCharge->Fill(ObsCharge);
                    hSimCharge->Fill(eLossTot);
                    hSim1Charge->Fill(eLoss1);
                    hSim2Charge->Fill(eLoss2);
                    hSim3Charge->Fill(eLoss3);
                    hSim4Charge->Fill(eLoss4);
                    hUnfoldCharge->Fill(UnfoldCharge);
                    hRecCharge->Fill(RecCharge);
                    hRecSelCharge->Fill(RecSelCharge);
                    hRelDiffCrossTalkWeighted->Fill(RelDiffCrossTalkWeighted);
                    hdEdx->Fill(GCCharge->at(c)/GCPath->at(c));
                    hRelDiffObs->Fill(RelDiffObs);
                    hRelDiffUnfold->Fill(RelDiffUnfold);
                    hRelDiffRec->Fill(RelDiffRec);
                    hRelDiffRecSel->Fill(RelDiffRecSel);
                    hRelDiffWeight->Fill(RelDiffWeight);

                }

                zlocaldist2=0;
                zlocaldist3=0;
                zlocaldist4=0;
                localdist2d2=0;
                localdist2d3=0;
                localdist2d4=0;
                localdist2=0;
                localdist3=0;
                localdist4=0;
                eLoss1=0.;
                eLoss2=0.;
                eLoss3=0.;
                eLoss4=0.;
                eLossTot=0.;
                thick=0.;
                ntracks=0.;
                wheel=0.;
                subdet=0.;
                q.clear();
                QUnfold.clear();
                QUnfold_P.clear();
                QUnfold_M.clear();
                QUnfoldProb.clear();
                Cluster.clear();
                Cluster13.clear();
                UnderCluster.clear();
                ExpCW=0.;
                UnfoldCW=0;
                Unfold13CW=0;
                RecCW=0.;
                RecSelCW=0.;
                RecSelCWUncert_P=0.;
                RecSelCWUncert_M=0.;
                RecSelProbCW=0.;
                RelDiff=0.;
                RelDiff1=0.;
                RelDiff12=0.;
                RelDiff2=0.;
                RelDiff23=0.;
                RelDiff3=0.;
                RelDiffCrossTalk_P=0.;
                RelDiffCrossTalk_M=0.;
                RelDiffCrossTalkProb=0.;
                RelDiffCrossTalkWeighted=0.;

            }
        }
      }


      binmax = hMeanRelDiffObs->GetMaximumBin();
      binx = hMeanRelDiffObs->GetXaxis()->GetBinCenter(binmax);
      hRelDiffThreshObs->SetBinContent(th+0.5,binx);
      hWidthThreshObs->SetBinContent(th+0.5,hMeanRelDiffObs->GetRMS());

      binmax = hMeanRelDiffUnfold->GetMaximumBin();
      binx = hMeanRelDiffUnfold->GetXaxis()->GetBinCenter(binmax);
      hRelDiffThreshUnfold->SetBinContent(th+0.5,binx);
      hWidthThreshUnfold->SetBinContent(th+0.5,hMeanRelDiffUnfold->GetRMS());

      binmax = hMeanRelDiffRec->GetMaximumBin();
      binx = hMeanRelDiffRec->GetXaxis()->GetBinCenter(binmax);
      hRelDiffThreshRec->SetBinContent(th+0.5,binx);
      hWidthThreshRec->SetBinContent(th+0.5,hMeanRelDiffRec->GetRMS());

      binmax = hMeanRelDiffRecSel->GetMaximumBin();
      binx = hMeanRelDiffRecSel->GetXaxis()->GetBinCenter(binmax);
      hRelDiffThreshRecSel->SetBinContent(th+0.5,binx);
      hWidthThreshRecSel->SetBinContent(th+0.5,hMeanRelDiffRecSel->GetRMS());

      binmax = hMeanRelDiffWeight->GetMaximumBin();
      binx = hMeanRelDiffWeight->GetXaxis()->GetBinCenter(binmax);
      hRelDiffThreshWeight->SetBinContent(th+0.5,binx);
      hWidthThreshWeight->SetBinContent(th+0.5, hMeanRelDiffWeight->GetRMS());

      hMeanRelDiffObs->Reset();
      hMeanRelDiffUnfold->Reset();
      hMeanRelDiffRec->Reset();
      hMeanRelDiffRecSel->Reset();
      hMeanRelDiffWeight->Reset();

    }

    TFile ofile(outfilename.c_str(),"RECREATE");

    hdeltaT->Write();
    hdeltaT13->Write();
    hdeltaT14->Write();

    hdeltaTdet->Write();
    hdeltaT13det->Write();
    hdeltaT14det->Write();

    hdeltaV->Write();

    hnTracks->Write();
    hsubdetratio->Write();
    hsubdet->Write();
    hsubdet1->Write();
    hsubdet2p->Write();

    hsimhits->Write();
    hdistOverThick2->Write();
    hdistOverThick3->Write();
    hdistOverThick4->Write();

    hzdistOverThick2->Write();
    hzdistOverThick3->Write();
    hzdistOverThick4->Write();

    h2DdistOverThick2->Write();
    h2DdistOverThick3->Write();
    h2DdistOverThick4->Write();

    hlocaldist2->Write();
    hlocaldist3->Write();
    hlocaldist4->Write();

    h2Dlocaldist2->Write();
    h2Dlocaldist3->Write();
    h2Dlocaldist4->Write();

    hzlocaldist2->Write();
    hzlocaldist3->Write();
    hzlocaldist4->Write();

    hwheelratio->Write();
    hwheel->Write();
    hwheel1->Write();
    hwheel2p->Write();

    hObsCW->Write();
    hExpCW->Write();
    hUnfoldCW->Write();
    hRecCW->Write();
    hRecSelCW->Write();
    hRecSelCWUncert_P->Write();
    hRecSelCWUncert_M->Write();
    hRelDiff->Write();
    hRelDiff1->Write();
    hRelDiff12->Write();
    hRelDiff2->Write();
    hRelDiff23->Write();
    hRelDiff3->Write();
    hRelDiffCrossTalk_P->Write();
    hRelDiffCrossTalk_M->Write();
    hWeightedWithUncertaintiesCrossTalkCW->Write();
    hRecSelProbCW->Write();
    hRelDiffCrossTalkProb->Write();
    hRelDiffCrossTalkWeighted->Write();
    hSimCharge->Write();
    hSim1Charge->Write();
    hSim2Charge->Write();
    hSim3Charge->Write();
    hSim4Charge->Write();
    hObsCharge->Write();
    hUnfoldCharge->Write();
    hRecCharge->Write();
    hRecSelCharge->Write();
    hdEdx->Write();

    hRelDiffObs->Write();
    hRelDiffUnfold->Write();
    hRelDiffRec->Write();
    hRelDiffRecSel->Write();
    hRelDiffWeight->Write();

    hRelDiffThreshObs->Write();
    hRelDiffThreshUnfold->Write();
    hRelDiffThreshRec->Write();
    hRelDiffThreshRecSel->Write();
    hRelDiffThreshWeight->Write();

    hWidthThreshObs->Write();
    hWidthThreshUnfold->Write();
    hWidthThreshRec->Write();
    hWidthThreshRecSel->Write();
    hWidthThreshWeight->Write();

    pObsCluster->Write();
    pUnfoldCluster->Write();


    hExpVsRec->GetXaxis()->SetTitle("Expected Cluster Width");
    hExpVsRec->GetYaxis()->SetTitle("Reconstructed Cluster Width");
    hExpVsRec->Write();
    c2->cd();
    gStyle->SetOptStat(0);
    hExpVsRec->Draw();
    hExpVsRec->SetDrawOption("COLZ");
    c2->SaveAs((plotDir+"ExpVsRec.pdf").c_str());
    c2->SaveAs((plotDir+"ExpVsRec.png").c_str());

    hExpVsObs->GetXaxis()->SetTitle("Expected Cluster Width");
    hExpVsObs->GetYaxis()->SetTitle("Observed Cluster Width");
    hExpVsObs->Write();
    c2->cd();
    gStyle->SetOptStat(0);
    hExpVsObs->Draw();
    hExpVsObs->SetDrawOption("COLZ");
    c2->SaveAs((plotDir+"ExpVsObs.pdf").c_str());
    c2->SaveAs((plotDir+"ExpVsObs.png").c_str());

    hExpVsRecSel->GetXaxis()->SetTitle("Expected Cluster Width");
    hExpVsRecSel->GetYaxis()->SetTitle("Reconstructed and Selected Cluster Width");
    hExpVsRecSel->Write();
    c2->cd();
    gStyle->SetOptStat(0);
    hExpVsRecSel->Draw();
    hExpVsRecSel->SetDrawOption("COLZ");
    c2->SaveAs((plotDir+"ExpVsRecSel.pdf").c_str());
    c2->SaveAs((plotDir+"ExpVsRecSel.png").c_str());

    hExpVsUnfold->GetXaxis()->SetTitle("Expected Cluster Width");
    hExpVsUnfold->GetYaxis()->SetTitle("Unfolded Cluster Width");
    hExpVsUnfold->Write();
    c2->cd();
    gStyle->SetOptStat(0);
    hExpVsUnfold->Draw();
    hExpVsUnfold->SetDrawOption("COLZ");
    c2->SaveAs((plotDir+"ExpVsUnfold.pdf").c_str());
    c2->SaveAs((plotDir+"ExpVsUnfold.png").c_str());

    hExpVsWeight->GetXaxis()->SetTitle("Expected Cluster Width");
    hExpVsWeight->GetYaxis()->SetTitle("Weighted Cluster Width");
    hExpVsWeight->Write();
    c2->cd();
    gStyle->SetOptStat(0);
    hExpVsWeight->Draw();
    hExpVsWeight->SetDrawOption("COLZ");
    c2->SaveAs((plotDir+"ExpVsWeight.pdf").c_str());
    c2->SaveAs((plotDir+"ExpVsWeight.png").c_str());

    hRelDiffCrossTalkWeighted->GetXaxis()->SetTitle("(Weighted - Normal)/Normal CrossTalk Cluster Width");
    c2->cd();
    gStyle->SetOptStat(0);
    hRelDiffCrossTalkWeighted->Draw();
    c2->SaveAs((plotDir+"RelDiffWeighted.pdf").c_str());
    c2->SaveAs((plotDir+"RelDiffWeighted.png").c_str());


    //--------- //ClusterWidth

    c1->cd();
    c1->SetLeftMargin(0.18);
    c1->SetBottomMargin(0.12);
    c1->SetLogy();
    hExpCW->GetXaxis()->SetTitle("Cluster Width");
    hExpCW->GetYaxis()->SetTitle("#Events");
//    hExpCW->GetYaxis()->SetRangeUser(100,40000);
    hExpCW->SetLineColor(38);
    hExpCW->SetLineWidth(3);
    hObsCW->SetLineColor(1);
    hObsCW->SetLineWidth(3);
    hUnfoldCW->SetLineColor(42);
    hUnfoldCW->SetLineWidth(3);
    hRecCW->SetLineColor(46);
    hRecCW->SetLineWidth(3);
    hRecSelCW->SetLineColor(8);
    hRecSelCW->SetLineWidth(3);
    hWeightedWithUncertaintiesCrossTalkCW->SetLineColor(40);
    hWeightedWithUncertaintiesCrossTalkCW->SetLineWidth(3);
    hRecSelProbCW->SetLineColor(28);
    hRecSelProbCW->SetLineWidth(3);

    hExpCW->GetXaxis()->SetLabelSize(0.045);
    hExpCW->GetYaxis()->SetLabelSize(0.045);
    hExpCW->GetXaxis()->SetTitleSize(0.045);
    hExpCW->GetYaxis()->SetTitleSize(0.045);

    hExpCW->Draw();
    hObsCW->Draw("SAME");
    hUnfoldCW->Draw("SAME");
    hRecCW->Draw("SAME");
//    hRecSelCW->Draw("SAME");
//    hWeightedWithUncertaintiesCrossTalkCW->Draw("SAME");
    // hRecSelProbCW->Draw("SAME");
    leg1->AddEntry(hObsCW,"Observed","l");
    leg1->AddEntry(hExpCW,"Expected","l");
    leg1->AddEntry(hUnfoldCW,"Unfolded","l");
    leg1->AddEntry(hRecCW,"Reconstructed","l");
//    leg1->AddEntry(hRecSelCW,"RecSel","l");
//    leg1->AddEntry(hWeightedWithUncertaintiesCrossTalkCW,"WeightedCrossTalk","l");
    // leg1->AddEntry(hRecSelProbCW,"RecSelProb","l");
    leg1->Draw("SAME");
    c1->Write();
    gStyle->SetOptStat(0);
    c1->SaveAs((plotDir+"ClusterWidth.pdf").c_str());
    c1->SaveAs((plotDir+"ClusterWidth.png").c_str());

    //--------- //Number of Tracks

    c5->cd();
    c5->SetLeftMargin(0.18);
    c5->SetBottomMargin(0.12);
    hnTracks->GetXaxis()->SetTitle("Number of Tracks per Cluster");
    hnTracks->SetLineColor(38);
    hnTracks->SetLineWidth(3);
    hnTracks->Draw();
    leg5->AddEntry(hnTracks,"Observed","l");
    leg5->Draw("SAME");
    c5->Write();
    gStyle->SetOptStat(0);
    c5->SaveAs((plotDir+"nTracks.pdf").c_str());
    c5->SaveAs((plotDir+"nTracks.png").c_str());

    //---------

    //--------- //Layer (subdetid)

    c6->cd();
    c6->SetLeftMargin(0.18);
    c6->SetBottomMargin(0.12);
    hsubdetratio->Divide(hsubdet);

    hsubdetratio->GetYaxis()->SetTitle("#Events (SimHits > 1) / #Events");
    hsubdetratio->GetXaxis()->SetTitle("Layer");
    hsubdetratio->SetLineColor(38);
    hsubdetratio->SetLineWidth(3);

    hsubdetratio->GetXaxis()->SetLabelSize(0.06);
    hsubdetratio->GetYaxis()->SetLabelSize(0.045);
    hsubdetratio->GetXaxis()->SetTitleSize(0.045);
    hsubdetratio->GetYaxis()->SetTitleSize(0.045);

    hsubdetratio->GetXaxis()->SetBinLabel( 1, "TIB" );
    hsubdetratio->GetXaxis()->SetBinLabel( 2, "TID" );
    hsubdetratio->GetXaxis()->SetBinLabel( 3, "TOB" );
    hsubdetratio->GetXaxis()->SetBinLabel( 4, "TEC" );

    hsubdetratio->Draw();
    leg6->AddEntry(hsubdetratio,"SimHits > 1 / total","l");
//    leg6->Draw("SAME");
    c6->Write();
    gStyle->SetOptStat(0);
    c6->SaveAs((plotDir+"layer_ratio.pdf").c_str());
    c6->SaveAs((plotDir+"layer_ratio.png").c_str());

    //---------

    //--------- //Layer (subdetid)

    c7->cd();
    c7->SetLeftMargin(0.18);
    c7->SetBottomMargin(0.12);

    hsubdet1->GetXaxis()->SetTitle("Layer");
    hsubdet1->GetYaxis()->SetTitle("#Events (SimHits == 1)");
    hsubdet1->SetLineColor(38);
    hsubdet1->SetLineWidth(3);

    hsubdet1->GetXaxis()->SetLabelSize(0.06);
    hsubdet1->GetYaxis()->SetLabelSize(0.045);
    hsubdet1->GetXaxis()->SetTitleSize(0.045);
    hsubdet1->GetYaxis()->SetTitleSize(0.045);

    hsubdet1->GetXaxis()->SetBinLabel( 1, "TIB" );
    hsubdet1->GetXaxis()->SetBinLabel( 2, "TID" );
    hsubdet1->GetXaxis()->SetBinLabel( 3, "TOB" );
    hsubdet1->GetXaxis()->SetBinLabel( 4, "TEC" );

    hsubdet1->Draw();
    leg7->AddEntry(hsubdet1,"SimHits = 1","l");
//    leg7->Draw("SAME");
    c7->Write();
    gStyle->SetOptStat(0);
    c7->SaveAs((plotDir+"layer_simhits1.pdf").c_str());
    c7->SaveAs((plotDir+"layer_simhits1.png").c_str());

    //---------

    //--------- //Layer (subdetid)

    c8->cd();
    c8->SetLeftMargin(0.18);
    c8->SetBottomMargin(0.12);

    hsubdet2p->GetYaxis()->SetTitle("#Events (SimHits > 1)");
    hsubdet2p->GetYaxis()->SetRangeUser(0,20000);
    hsubdet2p->GetXaxis()->SetTitle("Layer");
    hsubdet2p->SetLineColor(38);
    hsubdet2p->SetLineWidth(3);

    hsubdet2p->GetXaxis()->SetLabelSize(0.06);
    hsubdet2p->GetYaxis()->SetLabelSize(0.045);
    hsubdet2p->GetXaxis()->SetTitleSize(0.045);
    hsubdet2p->GetYaxis()->SetTitleSize(0.045);

    hsubdet2p->GetXaxis()->SetBinLabel( 1, "TIB" );
    hsubdet2p->GetXaxis()->SetBinLabel( 2, "TID" );
    hsubdet2p->GetXaxis()->SetBinLabel( 3, "TOB" );
    hsubdet2p->GetXaxis()->SetBinLabel( 4, "TEC" );

    hsubdet2p->Draw();
    leg8->AddEntry(hsubdet2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c8->Write();
    gStyle->SetOptStat(0);
    c8->SaveAs((plotDir+"layer_simhits2p.pdf").c_str());
    c8->SaveAs((plotDir+"layer_simhits2p.png").c_str());

    //---------



    //--------- //Layer (subdetid)

    c21->cd();
    c21->SetLeftMargin(0.18);
    c21->SetBottomMargin(0.12);

    hsubdet1->GetXaxis()->SetTitle("Layer");
    hsubdet1->GetYaxis()->SetTitle("Normalized #Events (SimHits == 1)");
    hsubdet1->SetLineColor(38);
    hsubdet1->SetLineWidth(3);

    hsubdet1->GetXaxis()->SetLabelSize(0.06);
    hsubdet1->GetYaxis()->SetLabelSize(0.045);
    hsubdet1->GetXaxis()->SetTitleSize(0.045);
    hsubdet1->GetYaxis()->SetTitleSize(0.045);

    hsubdet1->GetXaxis()->SetBinLabel( 1, "TIB" );
    hsubdet1->GetXaxis()->SetBinLabel( 2, "TID" );
    hsubdet1->GetXaxis()->SetBinLabel( 3, "TOB" );
    hsubdet1->GetXaxis()->SetBinLabel( 4, "TEC" );

    hsubdet1->Scale(1./hsubdet1->Integral());
    hsubdet1->Draw();
    leg7->AddEntry(hsubdet1,"SimHits = 1","l");
//    leg7->Draw("SAME");
    c21->Write();
    gStyle->SetOptStat(0);
    c21->SaveAs((plotDir+"layer_simhitsNorm1.pdf").c_str());
    c21->SaveAs((plotDir+"layer_simhitsNorm1.png").c_str());

    //---------

    //--------- //Layer (subdetid)

    c22->cd();
    c22->SetLeftMargin(0.18);
    c22->SetBottomMargin(0.12);

    hsubdet2p->GetYaxis()->SetTitle("Normalized #Events (SimHits > 1)");
    hsubdet2p->GetXaxis()->SetTitle("Layer");
    hsubdet2p->SetLineColor(38);
    hsubdet2p->SetLineWidth(3);

    hsubdet2p->GetXaxis()->SetLabelSize(0.06);
    hsubdet2p->GetYaxis()->SetLabelSize(0.045);
    hsubdet2p->GetXaxis()->SetTitleSize(0.045);
    hsubdet2p->GetYaxis()->SetTitleSize(0.045);

    hsubdet2p->GetXaxis()->SetBinLabel( 1, "TIB" );
    hsubdet2p->GetXaxis()->SetBinLabel( 2, "TID" );
    hsubdet2p->GetXaxis()->SetBinLabel( 3, "TOB" );
    hsubdet2p->GetXaxis()->SetBinLabel( 4, "TEC" );

    hsubdet2p->Scale(1./hsubdet2p->Integral());
    hsubdet2p->Draw();
    leg8->AddEntry(hsubdet2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c22->Write();
    gStyle->SetOptStat(0);
    c22->SaveAs((plotDir+"layer_simhitsNorm2p.pdf").c_str());
    c22->SaveAs((plotDir+"layer_simhitsNorm2p.png").c_str());

    //---------




    //--------- //time difference simhits

    c23->cd();
    c23->SetLeftMargin(0.15);
    c23->SetBottomMargin(0.15);

    hdeltaT->GetYaxis()->SetTitle("#Events");
    hdeltaT->GetXaxis()->SetTitle("#Delta t (simHit_{1}, simHit_{2})");
    hdeltaT->SetLineColor(38);
    hdeltaT->SetLineWidth(3);

    hdeltaT->GetXaxis()->SetLabelSize(0.045);
    hdeltaT->GetYaxis()->SetLabelSize(0.045);
    hdeltaT->GetXaxis()->SetTitleSize(0.045);
    hdeltaT->GetYaxis()->SetTitleSize(0.045);

    hdeltaT->Draw();
    c23->Write();
    gStyle->SetOptStat(0);
    c23->SaveAs((plotDir+"deltaT21.pdf").c_str());
    c23->SaveAs((plotDir+"deltaT21.png").c_str());

    hdeltaTdet->GetXaxis()->SetLabelSize(0.045);
    hdeltaTdet->GetYaxis()->SetLabelSize(0.045);
    hdeltaTdet->GetXaxis()->SetTitleSize(0.045);
    hdeltaTdet->GetYaxis()->SetTitleSize(0.045);

    hdeltaTdet->Draw();
    c23->SaveAs((plotDir+"deltaT21det.pdf").c_str());
    c23->SaveAs((plotDir+"deltaT21det.png").c_str());


    c25->cd();
    c25->SetLeftMargin(0.15);
    c25->SetBottomMargin(0.15);

    hdeltaT13->GetYaxis()->SetTitle("#Events");
    hdeltaT13->GetXaxis()->SetTitle("#Delta t (simHit_{1}, simHit_{3})");
    hdeltaT13->SetLineColor(38);
    hdeltaT13->SetLineWidth(3);

    hdeltaT13->GetXaxis()->SetLabelSize(0.045);
    hdeltaT13->GetYaxis()->SetLabelSize(0.045);
    hdeltaT13->GetXaxis()->SetTitleSize(0.045);
    hdeltaT13->GetYaxis()->SetTitleSize(0.045);

    hdeltaT13->Draw();
    c25->Write();
    gStyle->SetOptStat(0);
    c25->SaveAs((plotDir+"deltaT31.pdf").c_str());
    c25->SaveAs((plotDir+"deltaT31.png").c_str());

    hdeltaT13det->GetXaxis()->SetLabelSize(0.045);
    hdeltaT13det->GetYaxis()->SetLabelSize(0.045);
    hdeltaT13det->GetXaxis()->SetTitleSize(0.045);
    hdeltaT13det->GetYaxis()->SetTitleSize(0.045);

    hdeltaT13det->Draw();
    c25->SaveAs((plotDir+"deltaT31det.pdf").c_str());
    c25->SaveAs((plotDir+"deltaT31det.png").c_str());


    c26->cd();
    c26->SetLeftMargin(0.15);
    c26->SetBottomMargin(0.15);

    hdeltaT14->GetYaxis()->SetTitle("#Events");
    hdeltaT14->GetXaxis()->SetTitle("#Delta t (simHit_{1}, simHit_{4})");
    hdeltaT14->SetLineColor(38);
    hdeltaT14->SetLineWidth(3);

    hdeltaT14->GetXaxis()->SetLabelSize(0.045);
    hdeltaT14->GetYaxis()->SetLabelSize(0.045);
    hdeltaT14->GetXaxis()->SetTitleSize(0.045);
    hdeltaT14->GetYaxis()->SetTitleSize(0.045);

    hdeltaT14->Draw();
    c26->Write();
    gStyle->SetOptStat(0);
    c26->SaveAs((plotDir+"deltaT41.pdf").c_str());
    c26->SaveAs((plotDir+"deltaT41.png").c_str());

    hdeltaT14det->GetXaxis()->SetLabelSize(0.045);
    hdeltaT14det->GetYaxis()->SetLabelSize(0.045);
    hdeltaT14det->GetXaxis()->SetTitleSize(0.045);
    hdeltaT14det->GetYaxis()->SetTitleSize(0.045);

    hdeltaT14det->Draw();
    c26->SaveAs((plotDir+"deltaT41det.pdf").c_str());
    c26->SaveAs((plotDir+"deltaT41det.png").c_str());

    //---------



    //--------- //deltav simhits

    c24->cd();
    c24->SetLeftMargin(0.15);
    c24->SetBottomMargin(0.15);

    hdeltaV->GetYaxis()->SetTitle("#Events");
    hdeltaV->GetXaxis()->SetTitle("Sensor Thickness / #Delta t (simHit_{1}, simHit_{2})");
    hdeltaV->SetLineColor(38);
    hdeltaV->SetLineWidth(3);

    hdeltaV->GetXaxis()->SetLabelSize(0.045);
    hdeltaV->GetYaxis()->SetLabelSize(0.045);
    hdeltaV->GetXaxis()->SetTitleSize(0.045);
    hdeltaV->GetYaxis()->SetTitleSize(0.045);

    hdeltaV->Draw();
    c24->Write();
    gStyle->SetOptStat(0);
    c24->SaveAs((plotDir+"deltaV21.pdf").c_str());
    c24->SaveAs((plotDir+"deltaV21.png").c_str());

    //---------




    //--------- //Layer (wheel)

    c10->cd();
    c10->SetLeftMargin(0.18);
    c10->SetBottomMargin(0.12);
    hwheelratio->Divide(hwheel);

    hwheelratio->GetYaxis()->SetTitle("#Events (SimHits > 1) / #Events");
    hwheelratio->GetXaxis()->SetTitle("Layerwheel");
    hwheelratio->SetLineColor(38);
    hwheelratio->SetLineWidth(3);

    hwheelratio->GetXaxis()->SetLabelSize(0.06);
    hwheelratio->GetYaxis()->SetLabelSize(0.045);
    hwheelratio->GetXaxis()->SetTitleSize(0.045);
    hwheelratio->GetYaxis()->SetTitleSize(0.045);

    for (unsigned int w=1;w<10;w++){
        std::string s = std::to_string(w);
        hwheelratio->GetXaxis()->SetBinLabel( w+0.5, s.c_str() );
    }

    hwheelratio->Draw();
    leg6->AddEntry(hwheelratio,"SimHits > 1 / total","l");
//    leg6->Draw("SAME");
    c10->Write();
    gStyle->SetOptStat(0);
    c10->SaveAs((plotDir+"wheel_ratio.pdf").c_str());
    c10->SaveAs((plotDir+"wheel_ratio.png").c_str());

    //---------

    //--------- //Layer (wheelid)

    c11->cd();
    c11->SetLeftMargin(0.18);
    c11->SetBottomMargin(0.12);

    hwheel1->GetXaxis()->SetTitle("Layerwheel");
    hwheel1->GetYaxis()->SetTitle("#Events (SimHits == 1)");
    hwheel1->SetLineColor(38);
    hwheel1->SetLineWidth(3);

    hwheel1->GetXaxis()->SetLabelSize(0.06);
    hwheel1->GetYaxis()->SetLabelSize(0.045);
    hwheel1->GetXaxis()->SetTitleSize(0.045);
    hwheel1->GetYaxis()->SetTitleSize(0.045);

    for (unsigned int w=1;w<10;w++){
        std::string s = std::to_string(w);
        hwheel1->GetXaxis()->SetBinLabel( w+0.5, s.c_str() );
    }

    hwheel1->Draw();
    leg7->AddEntry(hwheel1,"SimHits = 1","l");
//    leg7->Draw("SAME");
    c11->Write();
    gStyle->SetOptStat(0);
    c11->SaveAs((plotDir+"wheel_simhits1.pdf").c_str());
    c11->SaveAs((plotDir+"wheel_simhits1.png").c_str());

    //---------

    //--------- //Layer (wheelid)

    c12->cd();
    c12->SetLeftMargin(0.18);
    c12->SetBottomMargin(0.12);

    hwheel2p->GetYaxis()->SetTitle("#Events (SimHits > 1)");
    hwheel2p->GetXaxis()->SetTitle("Layerwheel");
    hwheel2p->SetLineColor(38);
    hwheel2p->SetLineWidth(3);

    hwheel2p->GetXaxis()->SetLabelSize(0.06);
    hwheel2p->GetYaxis()->SetLabelSize(0.045);
    hwheel2p->GetXaxis()->SetTitleSize(0.045);
    hwheel2p->GetYaxis()->SetTitleSize(0.045);

    for (unsigned int w=1;w<10;w++){
        std::string s = std::to_string(w);
        hwheel2p->GetXaxis()->SetBinLabel( w+0.5, s.c_str() );
    }

    hwheel2p->Draw();
    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c12->Write();
    gStyle->SetOptStat(0);
    c12->SaveAs((plotDir+"wheel_simhits2p.pdf").c_str());
    c12->SaveAs((plotDir+"wheel_simhits2p.png").c_str());

    //---------


    //--------- //simhits

    c13->cd();
    c13->SetLeftMargin(0.18);
    c13->SetBottomMargin(0.12);
    c13->SetLogy();

    hsimhits->GetYaxis()->SetTitle("#Events");
    hsimhits->GetXaxis()->SetTitle("Simulated Hits");
    hsimhits->SetLineColor(38);
    hsimhits->SetLineWidth(3);

    hsimhits->GetXaxis()->SetLabelSize(0.06);
    hsimhits->GetYaxis()->SetLabelSize(0.045);
    hsimhits->GetXaxis()->SetTitleSize(0.045);
    hsimhits->GetYaxis()->SetTitleSize(0.045);

    hsimhits->GetXaxis()->SetBinLabel( 1, "1 hit" );
    hsimhits->GetXaxis()->SetBinLabel( 2, "2 hits" );
    hsimhits->GetXaxis()->SetBinLabel( 3, "3 hits" );

    hsimhits->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c13->Write();
    gStyle->SetOptStat(0);
    c13->SaveAs((plotDir+"simhits.pdf").c_str());
    c13->SaveAs((plotDir+"simhits.png").c_str());

    //---------


    //--------- //simhits

    c14->cd();
    c14->SetLeftMargin(0.15);
    c14->SetBottomMargin(0.15);

    hlocaldist2->GetYaxis()->SetTitle("#Events");
    hlocaldist2->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{2})");
    hlocaldist2->SetLineColor(38);
    hlocaldist2->SetLineWidth(3);

    hlocaldist2->GetXaxis()->SetLabelSize(0.045);
    hlocaldist2->GetYaxis()->SetLabelSize(0.045);
    hlocaldist2->GetXaxis()->SetTitleSize(0.045);
    hlocaldist2->GetYaxis()->SetTitleSize(0.045);

    hlocaldist2->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c14->Write();
    gStyle->SetOptStat(0);
    c14->SaveAs((plotDir+"simhits_dist12.pdf").c_str());
    c14->SaveAs((plotDir+"simhits_dist12.png").c_str());



    c16->cd();
    c16->SetLeftMargin(0.15);
    c16->SetBottomMargin(0.15);

    hlocaldist3->GetYaxis()->SetTitle("#Events");
    hlocaldist3->GetXaxis()->SetTitle("Rel. Dist. (simHit_{1}, simHit_{3})");
    hlocaldist3->SetLineColor(38);
    hlocaldist3->SetLineWidth(3);

    hlocaldist3->GetXaxis()->SetLabelSize(0.045);
    hlocaldist3->GetYaxis()->SetLabelSize(0.045);
    hlocaldist3->GetXaxis()->SetTitleSize(0.045);
    hlocaldist3->GetYaxis()->SetTitleSize(0.045);

    hlocaldist3->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c16->Write();
    gStyle->SetOptStat(0);
    c16->SaveAs((plotDir+"simhits_dist13.pdf").c_str());
    c16->SaveAs((plotDir+"simhits_dist13.png").c_str());



    c17->cd();
    c17->SetLeftMargin(0.15);
    c17->SetBottomMargin(0.15);

    hlocaldist4->GetYaxis()->SetTitle("#Events");
    hlocaldist4->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{4})");
    hlocaldist4->SetLineColor(38);
    hlocaldist4->SetLineWidth(3);

    hlocaldist4->GetXaxis()->SetLabelSize(0.045);
    hlocaldist4->GetYaxis()->SetLabelSize(0.045);
    hlocaldist4->GetXaxis()->SetTitleSize(0.045);
    hlocaldist4->GetYaxis()->SetTitleSize(0.045);

    hlocaldist4->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c17->Write();
    gStyle->SetOptStat(0);
    c17->SaveAs((plotDir+"simhits_dist14.pdf").c_str());
    c17->SaveAs((plotDir+"simhits_dist14.png").c_str());


    c18->cd();
    c18->SetLeftMargin(0.15);
    c18->SetBottomMargin(0.15);

    hdistOverThick2->GetYaxis()->SetTitle("#Events");
    hdistOverThick2->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{2}) / Sensor thickness");
    hdistOverThick2->SetLineColor(38);
    hdistOverThick2->SetLineWidth(3);

    hdistOverThick2->GetXaxis()->SetLabelSize(0.045);
    hdistOverThick2->GetYaxis()->SetLabelSize(0.045);
    hdistOverThick2->GetXaxis()->SetTitleSize(0.045);
    hdistOverThick2->GetYaxis()->SetTitleSize(0.045);

    hdistOverThick2->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c18->Write();
    gStyle->SetOptStat(0);
    c18->SaveAs((plotDir+"simhits_distOverThick12.pdf").c_str());
    c18->SaveAs((plotDir+"simhits_distOverThick12.png").c_str());


    c19->cd();
    c19->SetLeftMargin(0.15);
    c19->SetBottomMargin(0.15);

    hdistOverThick3->GetYaxis()->SetTitle("#Events");
    hdistOverThick3->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{3}) / Sensor thickness");
    hdistOverThick3->SetLineColor(38);
    hdistOverThick3->SetLineWidth(3);

    hdistOverThick3->GetXaxis()->SetLabelSize(0.045);
    hdistOverThick3->GetYaxis()->SetLabelSize(0.045);
    hdistOverThick3->GetXaxis()->SetTitleSize(0.045);
    hdistOverThick3->GetYaxis()->SetTitleSize(0.045);

    hdistOverThick3->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c19->Write();
    gStyle->SetOptStat(0);
    c19->SaveAs((plotDir+"simhits_distOverThick13.pdf").c_str());
    c19->SaveAs((plotDir+"simhits_distOverThick13.png").c_str());


    c20->cd();
    c20->SetLeftMargin(0.15);
    c20->SetBottomMargin(0.15);

    hdistOverThick4->GetYaxis()->SetTitle("#Events");
    hdistOverThick4->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{4}) / Sensor thickness");
    hdistOverThick4->SetLineColor(38);
    hdistOverThick4->SetLineWidth(3);

    hdistOverThick4->GetXaxis()->SetLabelSize(0.045);
    hdistOverThick4->GetYaxis()->SetLabelSize(0.045);
    hdistOverThick4->GetXaxis()->SetTitleSize(0.045);
    hdistOverThick4->GetYaxis()->SetTitleSize(0.045);

    hdistOverThick4->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c20->Write();
    gStyle->SetOptStat(0);
    c20->SaveAs((plotDir+"simhits_distOverThick14.pdf").c_str());
    c20->SaveAs((plotDir+"simhits_distOverThick14.png").c_str());

    //---------



    c27->cd();
    c27->SetLeftMargin(0.15);
    c27->SetBottomMargin(0.15);

    h2Dlocaldist2->GetYaxis()->SetTitle("#Events");
    h2Dlocaldist2->GetXaxis()->SetTitle("2D Rel. dist. (simHit_{1}, simHit_{2})");
    h2Dlocaldist2->SetLineColor(38);
    h2Dlocaldist2->SetLineWidth(3);

    h2Dlocaldist2->GetXaxis()->SetLabelSize(0.045);
    h2Dlocaldist2->GetYaxis()->SetLabelSize(0.045);
    h2Dlocaldist2->GetXaxis()->SetTitleSize(0.045);
    h2Dlocaldist2->GetYaxis()->SetTitleSize(0.045);

    h2Dlocaldist2->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c27->Write();
    gStyle->SetOptStat(0);
    c27->SaveAs((plotDir+"simhits_2Ddist12.pdf").c_str());
    c27->SaveAs((plotDir+"simhits_2Ddist12.png").c_str());



    c28->cd();
    c28->SetLeftMargin(0.15);
    c28->SetBottomMargin(0.15);

    h2Dlocaldist3->GetYaxis()->SetTitle("#Events");
    h2Dlocaldist3->GetXaxis()->SetTitle("2D Rel. Dist. (simHit_{1}, simHit_{3})");
    h2Dlocaldist3->SetLineColor(38);
    h2Dlocaldist3->SetLineWidth(3);

    h2Dlocaldist3->GetXaxis()->SetLabelSize(0.045);
    h2Dlocaldist3->GetYaxis()->SetLabelSize(0.045);
    h2Dlocaldist3->GetXaxis()->SetTitleSize(0.045);
    h2Dlocaldist3->GetYaxis()->SetTitleSize(0.045);

    h2Dlocaldist3->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c28->Write();
    gStyle->SetOptStat(0);
    c28->SaveAs((plotDir+"simhits_2Ddist13.pdf").c_str());
    c28->SaveAs((plotDir+"simhits_2Ddist13.png").c_str());



    c29->cd();
    c29->SetLeftMargin(0.15);
    c29->SetBottomMargin(0.15);

    h2Dlocaldist4->GetYaxis()->SetTitle("#Events");
    h2Dlocaldist4->GetXaxis()->SetTitle("2D Rel. dist. (simHit_{1}, simHit_{4})");
    h2Dlocaldist4->SetLineColor(38);
    h2Dlocaldist4->SetLineWidth(3);

    h2Dlocaldist4->GetXaxis()->SetLabelSize(0.045);
    h2Dlocaldist4->GetYaxis()->SetLabelSize(0.045);
    h2Dlocaldist4->GetXaxis()->SetTitleSize(0.045);
    h2Dlocaldist4->GetYaxis()->SetTitleSize(0.045);

    h2Dlocaldist4->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c29->Write();
    gStyle->SetOptStat(0);
    c29->SaveAs((plotDir+"simhits_2Ddist14.pdf").c_str());
    c29->SaveAs((plotDir+"simhits_2Ddist14.png").c_str());


    c30->cd();
    c30->SetLeftMargin(0.15);
    c30->SetBottomMargin(0.15);

    h2DdistOverThick2->GetYaxis()->SetTitle("#Events");
    h2DdistOverThick2->GetXaxis()->SetTitle("2D Rel. dist. (simHit_{1}, simHit_{2}) / Sensor thickness");
    h2DdistOverThick2->SetLineColor(38);
    h2DdistOverThick2->SetLineWidth(3);

    h2DdistOverThick2->GetXaxis()->SetLabelSize(0.045);
    h2DdistOverThick2->GetYaxis()->SetLabelSize(0.045);
    h2DdistOverThick2->GetXaxis()->SetTitleSize(0.045);
    h2DdistOverThick2->GetYaxis()->SetTitleSize(0.045);

    h2DdistOverThick2->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c30->Write();
    gStyle->SetOptStat(0);
    c30->SaveAs((plotDir+"simhits_2DdistOverThick12.pdf").c_str());
    c30->SaveAs((plotDir+"simhits_2DdistOverThick12.png").c_str());


    c31->cd();
    c31->SetLeftMargin(0.15);
    c31->SetBottomMargin(0.15);

    h2DdistOverThick3->GetYaxis()->SetTitle("#Events");
    h2DdistOverThick3->GetXaxis()->SetTitle("2D Rel. dist. (simHit_{1}, simHit_{3}) / Sensor thickness");
    h2DdistOverThick3->SetLineColor(38);
    h2DdistOverThick3->SetLineWidth(3);

    h2DdistOverThick3->GetXaxis()->SetLabelSize(0.045);
    h2DdistOverThick3->GetYaxis()->SetLabelSize(0.045);
    h2DdistOverThick3->GetXaxis()->SetTitleSize(0.045);
    h2DdistOverThick3->GetYaxis()->SetTitleSize(0.045);

    h2DdistOverThick3->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c31->Write();
    gStyle->SetOptStat(0);
    c31->SaveAs((plotDir+"simhits_2DdistOverThick13.pdf").c_str());
    c31->SaveAs((plotDir+"simhits_2DdistOverThick13.png").c_str());


    c32->cd();
    c32->SetLeftMargin(0.15);
    c32->SetBottomMargin(0.15);

    h2DdistOverThick4->GetYaxis()->SetTitle("#Events");
    h2DdistOverThick4->GetXaxis()->SetTitle("2D Rel. dist. (simHit_{1}, simHit_{4}) / Sensor thickness");
    h2DdistOverThick4->SetLineColor(38);
    h2DdistOverThick4->SetLineWidth(3);

    h2DdistOverThick4->GetXaxis()->SetLabelSize(0.045);
    h2DdistOverThick4->GetYaxis()->SetLabelSize(0.045);
    h2DdistOverThick4->GetXaxis()->SetTitleSize(0.045);
    h2DdistOverThick4->GetYaxis()->SetTitleSize(0.045);

    h2DdistOverThick4->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c32->Write();
    gStyle->SetOptStat(0);
    c32->SaveAs((plotDir+"simhits_2DdistOverThick14.pdf").c_str());
    c32->SaveAs((plotDir+"simhits_2DdistOverThick14.png").c_str());

    //---------



    c33->cd();
    c33->SetLeftMargin(0.15);
    c33->SetBottomMargin(0.15);

    hzlocaldist2->GetYaxis()->SetTitle("#Events");
    hzlocaldist2->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{2}) (z)");
    hzlocaldist2->SetLineColor(38);
    hzlocaldist2->SetLineWidth(3);

    hzlocaldist2->GetXaxis()->SetLabelSize(0.045);
    hzlocaldist2->GetYaxis()->SetLabelSize(0.045);
    hzlocaldist2->GetXaxis()->SetTitleSize(0.045);
    hzlocaldist2->GetYaxis()->SetTitleSize(0.045);

    hzlocaldist2->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c33->Write();
    gStyle->SetOptStat(0);
    c33->SaveAs((plotDir+"simhits_zdist12.pdf").c_str());
    c33->SaveAs((plotDir+"simhits_zdist12.png").c_str());



    c34->cd();
    c34->SetLeftMargin(0.15);
    c34->SetBottomMargin(0.15);

    hzlocaldist3->GetYaxis()->SetTitle("#Events");
    hzlocaldist3->GetXaxis()->SetTitle("Rel. Dist. (simHit_{1}, simHit_{3}) (z)");
    hzlocaldist3->SetLineColor(38);
    hzlocaldist3->SetLineWidth(3);

    hzlocaldist3->GetXaxis()->SetLabelSize(0.045);
    hzlocaldist3->GetYaxis()->SetLabelSize(0.045);
    hzlocaldist3->GetXaxis()->SetTitleSize(0.045);
    hzlocaldist3->GetYaxis()->SetTitleSize(0.045);

    hzlocaldist3->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c34->Write();
    gStyle->SetOptStat(0);
    c34->SaveAs((plotDir+"simhits_zdist13.pdf").c_str());
    c34->SaveAs((plotDir+"simhits_zdist13.png").c_str());



    c35->cd();
    c35->SetLeftMargin(0.15);
    c35->SetBottomMargin(0.15);

    hzlocaldist4->GetYaxis()->SetTitle("#Events");
    hzlocaldist4->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{4}) (z)");
    hzlocaldist4->SetLineColor(38);
    hzlocaldist4->SetLineWidth(3);

    hzlocaldist4->GetXaxis()->SetLabelSize(0.045);
    hzlocaldist4->GetYaxis()->SetLabelSize(0.045);
    hzlocaldist4->GetXaxis()->SetTitleSize(0.045);
    hzlocaldist4->GetYaxis()->SetTitleSize(0.045);

    hzlocaldist4->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c35->Write();
    gStyle->SetOptStat(0);
    c35->SaveAs((plotDir+"simhits_zdist14.pdf").c_str());
    c35->SaveAs((plotDir+"simhits_zdist14.png").c_str());



    c36->cd();
    c36->SetLeftMargin(0.15);
    c36->SetBottomMargin(0.15);

    hzdistOverThick2->GetYaxis()->SetTitle("#Events");
    hzdistOverThick2->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{2}) (z) / Sensor thickness");
    hzdistOverThick2->SetLineColor(38);
    hzdistOverThick2->SetLineWidth(3);

    hzdistOverThick2->GetXaxis()->SetLabelSize(0.045);
    hzdistOverThick2->GetYaxis()->SetLabelSize(0.045);
    hzdistOverThick2->GetXaxis()->SetTitleSize(0.045);
    hzdistOverThick2->GetYaxis()->SetTitleSize(0.045);

    hzdistOverThick2->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c36->Write();
    gStyle->SetOptStat(0);
    c36->SaveAs((plotDir+"simhits_zdistOverThick12.pdf").c_str());
    c36->SaveAs((plotDir+"simhits_zdistOverThick12.png").c_str());


    c37->cd();
    c37->SetLeftMargin(0.15);
    c37->SetBottomMargin(0.15);

    hzdistOverThick3->GetYaxis()->SetTitle("#Events");
    hzdistOverThick3->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{3}) (z) / Sensor thickness");
    hzdistOverThick3->SetLineColor(38);
    hzdistOverThick3->SetLineWidth(3);

    hzdistOverThick3->GetXaxis()->SetLabelSize(0.045);
    hzdistOverThick3->GetYaxis()->SetLabelSize(0.045);
    hzdistOverThick3->GetXaxis()->SetTitleSize(0.045);
    hzdistOverThick3->GetYaxis()->SetTitleSize(0.045);

    hzdistOverThick3->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c37->Write();
    gStyle->SetOptStat(0);
    c37->SaveAs((plotDir+"simhits_zdistOverThick13.pdf").c_str());
    c37->SaveAs((plotDir+"simhits_zdistOverThick13.png").c_str());


    c38->cd();
    c38->SetLeftMargin(0.15);
    c38->SetBottomMargin(0.15);

    hzdistOverThick4->GetYaxis()->SetTitle("#Events");
    hzdistOverThick4->GetXaxis()->SetTitle("Rel. dist. (simHit_{1}, simHit_{4}) (z) / Sensor thickness");
    hzdistOverThick4->SetLineColor(38);
    hzdistOverThick4->SetLineWidth(3);

    hzdistOverThick4->GetXaxis()->SetLabelSize(0.045);
    hzdistOverThick4->GetYaxis()->SetLabelSize(0.045);
    hzdistOverThick4->GetXaxis()->SetTitleSize(0.045);
    hzdistOverThick4->GetYaxis()->SetTitleSize(0.045);

    hzdistOverThick4->Draw();
//    leg8->AddEntry(hwheel2p,"SimHits > 1","l");
//    leg8->Draw("SAME");
    c38->Write();
    gStyle->SetOptStat(0);
    c38->SaveAs((plotDir+"simhits_zdistOverThick14.pdf").c_str());
    c38->SaveAs((plotDir+"simhits_zdistOverThick14.png").c_str());

    //---------

    if(PrintThresholds){
        //--------- //RelDiff in thresh

        c9->cd();
        c9->SetLeftMargin(0.15);
        c9->SetBottomMargin(0.15);

        TLegend* legRelDiffThresh=new TLegend(0.55,0.65,0.9,0.9);
        c9->cd();
        hRelDiffThreshObs->GetYaxis()->SetTitle("MPV( (Rec-Exp)/Exp Cluster Width )");
        hRelDiffThreshObs->GetXaxis()->SetTitle("Threshold");
        hRelDiffThreshObs->SetLineColor(1);
        hRelDiffThreshObs->SetLineWidth(3);
        hRelDiffThreshUnfold->SetLineColor(42);
        hRelDiffThreshUnfold->SetLineWidth(3);
        hRelDiffThreshRec->SetLineColor(46);
        hRelDiffThreshRec->SetLineWidth(3);
        hRelDiffThreshRecSel->SetLineColor(8);
        hRelDiffThreshRecSel->SetLineWidth(3);
        hRelDiffThreshWeight->SetLineColor(40);
        hRelDiffThreshWeight->SetLineWidth(3);
        hRelDiffThreshObs->GetYaxis()->SetRangeUser(0,1);

        hRelDiffThreshObs->GetXaxis()->SetLabelSize(0.045);
        hRelDiffThreshObs->GetYaxis()->SetLabelSize(0.045);
        hRelDiffThreshObs->GetXaxis()->SetTitleSize(0.045);
        hRelDiffThreshObs->GetYaxis()->SetTitleSize(0.045);

        hRelDiffThreshObs->Draw();
        hRelDiffThreshUnfold->Draw("SAME");
        hRelDiffThreshRec->Draw("SAME");
//        hRelDiffThreshRecSel->Draw("SAME");
//        hRelDiffThreshWeight->Draw("SAME");
        legRelDiffThresh->AddEntry(hRelDiffThreshObs,"Observed","l");
        legRelDiffThresh->AddEntry(hRelDiffThreshUnfold,"Unfolded","l");
        legRelDiffThresh->AddEntry(hRelDiffThreshRec,"Reconstructed","l");
//        legRelDiffThresh->AddEntry(hRelDiffThreshRecSel,"RecSel","l");
//        legRelDiffThresh->AddEntry(hRelDiffThreshWeight,"Weight","l");
        legRelDiffThresh->Draw("SAME");
        c9->Write();
        c9->SaveAs((plotDir+"Thresh_RelDiff.pdf").c_str());
        c9->SaveAs((plotDir+"Thresh_RelDiff.png").c_str());
    //---------


        //--------- //RelDiff in thresh

        c15->cd();
        c15->SetLeftMargin(0.15);
        c15->SetBottomMargin(0.15);

        TLegend* legWidthThresh=new TLegend(0.55,0.65,0.9,0.9);
        c15->cd();
        hWidthThreshObs->GetYaxis()->SetTitle("Std. Dev.( (Rec-Exp)/Exp Cluster Width )");
        hWidthThreshObs->GetXaxis()->SetTitle("Threshold");
        hWidthThreshObs->SetLineColor(1);
        hWidthThreshObs->SetLineWidth(3);
        hWidthThreshUnfold->SetLineColor(42);
        hWidthThreshUnfold->SetLineWidth(3);
        hWidthThreshRec->SetLineColor(46);
        hWidthThreshRec->SetLineWidth(3);
        hWidthThreshRecSel->SetLineColor(8);
        hWidthThreshRecSel->SetLineWidth(3);
        hWidthThreshWeight->SetLineColor(40);
        hWidthThreshWeight->SetLineWidth(3);
        hWidthThreshObs->GetYaxis()->SetRangeUser(0,1);

        hWidthThreshObs->GetXaxis()->SetLabelSize(0.045);
        hWidthThreshObs->GetYaxis()->SetLabelSize(0.045);
        hWidthThreshObs->GetXaxis()->SetTitleSize(0.045);
        hWidthThreshObs->GetYaxis()->SetTitleSize(0.045);

        hWidthThreshObs->Draw();
        hWidthThreshUnfold->Draw("SAME");
        hWidthThreshRec->Draw("SAME");
//        hWidthThreshRecSel->Draw("SAME");
//        hWidthThreshWeight->Draw("SAME");
        legWidthThresh->AddEntry(hWidthThreshObs,"Observed","l");
        legWidthThresh->AddEntry(hWidthThreshUnfold,"Unfolded","l");
        legWidthThresh->AddEntry(hWidthThreshRec,"Reconstructed","l");
//        legWidthThresh->AddEntry(hWidthThreshRecSel,"RecSel","l");
//        legWidthThresh->AddEntry(hWidthThreshWeight,"Weight","l");
        legWidthThresh->Draw("SAME");
        c15->Write();
        c15->SaveAs((plotDir+"Thresh_Width.pdf").c_str());
        c15->SaveAs((plotDir+"Thresh_Width.png").c_str());
    //---------
    }

    //--------- //Charge
    c3->cd();
    c3->SetLeftMargin(0.18);
    c3->SetBottomMargin(0.12);

    hSimCharge->SetLineColor(8);
    hSimCharge->SetLineWidth(3);
    hObsCharge->SetLineColor(1);
    hObsCharge->SetLineWidth(3);
    hUnfoldCharge->SetLineColor(38);
    hUnfoldCharge->SetLineWidth(3);
    hRecCharge->SetLineColor(46);
    hRecCharge->SetLineWidth(3);
    hRecSelCharge->SetLineColor(8);
    hRecSelCharge->SetLineWidth(3);
    hSimCharge->GetYaxis()->SetTitle("#Events");
    hSimCharge->GetXaxis()->SetTitle("Charge");

    hSimCharge->GetXaxis()->SetLabelSize(0.045);
    hSimCharge->GetYaxis()->SetLabelSize(0.045);
    hSimCharge->GetXaxis()->SetTitleSize(0.045);
    hSimCharge->GetYaxis()->SetTitleSize(0.045);

    hSimCharge->Draw();
    hObsCharge->Draw("SAME");
    hUnfoldCharge->Draw("SAME");
    hRecCharge->Draw("SAME");
//    hRecSelCharge->Draw("SAME");
    leg4->AddEntry(hSimCharge,"Simulated","l");
    leg4->AddEntry(hObsCharge,"Observed","l");
    leg4->AddEntry(hUnfoldCharge,"Unfolded","l");
    leg4->AddEntry(hRecCharge,"Reconstructed","l");
//    leg4->AddEntry(hRecSelCharge,"RecSel","l");
    leg4->Draw("SAME");
    c3->Write();
    c3->SaveAs((plotDir+"Charge.pdf").c_str());
    c3->SaveAs((plotDir+"Charge.png").c_str());
    //---------

    //--------- //Charge
    c39->cd();
    c39->SetLeftMargin(0.18);
    c39->SetBottomMargin(0.12);
    hSimCharge->SetLineColor(2);
    hSimCharge->SetLineWidth(3);
    hSim1Charge->SetLineColor(1);
    hSim1Charge->SetLineWidth(3);
    hSim2Charge->SetLineColor(38);
    hSim2Charge->SetLineWidth(3);
    hSim3Charge->SetLineColor(46);
    hSim3Charge->SetLineWidth(3);
    hSim4Charge->SetLineColor(8);
    hSim4Charge->SetLineWidth(3);
    hSim1Charge->GetYaxis()->SetTitle("#Events");
    hSim1Charge->GetXaxis()->SetTitle("Charge");
//    hSimCharge->Draw();

    hSim1Charge->GetXaxis()->SetLabelSize(0.045);
    hSim1Charge->GetYaxis()->SetLabelSize(0.045);
    hSim1Charge->GetXaxis()->SetTitleSize(0.045);
    hSim1Charge->GetYaxis()->SetTitleSize(0.045);

    hSim1Charge->Draw();
    hSim2Charge->Draw("SAME");
    hSim3Charge->Draw("SAME");
    hSim4Charge->Draw("SAME");
//    leg10->AddEntry(hSimCharge,"Simulated (comb)","l");
    leg10->AddEntry(hSim1Charge,"Simulated (hit 1)","l");
    leg10->AddEntry(hSim2Charge,"Simulated (hit 2)","l");
    leg10->AddEntry(hSim3Charge,"Simulated (hit 3)","l");
    leg10->AddEntry(hSim4Charge,"Simulated (hit 4)","l");
    leg10->Draw("SAME");
    c39->Write();
    c39->SaveAs((plotDir+"SimCharge.pdf").c_str());
    c39->SaveAs((plotDir+"SimCharge.png").c_str());
    //---------

    //--------- //RelDiff

    TLegend* legRelDiff=new TLegend(0.55,0.65,0.9,0.9);
    c2->cd();
    c2->SetLeftMargin(0.18);
    c2->SetBottomMargin(0.12);
    hRelDiffObs->GetYaxis()->SetTitle("#Events");
    hRelDiffObs->GetXaxis()->SetTitle("(CW - ExpCW)/ExpCW");
    hRelDiffObs->GetYaxis()->SetRangeUser(0,hRelDiffRec->GetMaximum()+9);
    hRelDiffObs->SetLineColor(1);
    hRelDiffObs->SetLineWidth(3);
    hRelDiffUnfold->SetLineColor(42);
    hRelDiffUnfold->SetLineWidth(3);
    hRelDiffRec->SetLineColor(46);
    hRelDiffRec->SetLineWidth(3);
    hRelDiffRecSel->SetLineColor(8);
    hRelDiffRecSel->SetLineWidth(3);
    hRelDiffWeight->SetLineColor(40);
    hRelDiffWeight->SetLineWidth(3);

    hRelDiffObs->GetXaxis()->SetLabelSize(0.045);
    hRelDiffObs->GetYaxis()->SetLabelSize(0.045);
    hRelDiffObs->GetXaxis()->SetTitleSize(0.045);
    hRelDiffObs->GetYaxis()->SetTitleSize(0.045);

    hRelDiffObs->Draw();
    hRelDiffUnfold->Draw("SAME");
    hRelDiffRec->Draw("SAME");
//    hRelDiffRecSel->Draw("SAME");
//    hRelDiffWeight->Draw("SAME");
    legRelDiff->AddEntry(hRelDiffObs,"Observed","l");
    legRelDiff->AddEntry(hRelDiffUnfold,"Unfolded","l");
    legRelDiff->AddEntry(hRelDiffRec,"Reconstructed","l");
//    legRelDiff->AddEntry(hRelDiffRecSel,"RecSel","l");
//    legRelDiff->AddEntry(hRelDiffWeight,"Weight","l");
    legRelDiff->Draw("SAME");
    c2->Write();
    c2->SaveAs((plotDir+"RelDiff2.pdf").c_str());
    c2->SaveAs((plotDir+"RelDiff2.png").c_str());



    //---------


    ofile.Write();
    ofile.Close();

    return 0;
}
