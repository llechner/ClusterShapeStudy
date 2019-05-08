#include <iostream>
#include <vector>
#include <list>
#include <string>
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

#include "../interface/Estimator.h"
#include "../interface/ClusterShape.h"

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

const double clusterThreshold=13.; // based on observation, make that flexible?

using namespace std;

int main(){

	//------------------------------------ //Ouverture fichier+SetBranch

	// const char *path="/eos/cms/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR17_Aag/";
	//string filename="/eos/cms/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR17_Aag/calibTree_304777_25.root";
	string filename="/afs/cern.ch/work/l/llechner/public/CMSSW_9_4_4/src/cmssw/ClusterShape/data/merged_1.root";
	string outfilename="/afs/cern.ch/work/l/llechner/public/CMSSW_9_4_4/src/cmssw/ClusterShape/data/resutls.root";

	TFile* ifile=TFile::Open(filename.c_str());
	TTree* tree=(TTree*) ifile->Get("testTree/tree");

	std::vector<unsigned int>*clusterStripIdx	= 0;
	std::vector<double>* clusterAmplitudes		= 0;
	std::vector<unsigned int>*clusterIdx		= 0;
	std::vector<double>* clusterCharge		= 0;
	std::vector<unsigned int>*clusterWidth		= 0;
	std::vector<double>* clusterLocalTrackPhi	= 0;
	std::vector<double>* clusterLocalTrackTheta	= 0;
	std::vector<int>*clusterSubdetid		= 0;
	std::vector<int>*clusterLayerwheel		= 0;
	std::vector<unsigned int>*clusterDetid		= 0;
	std::vector<double>* clusterLocalTrackX		= 0;
	std::vector<double>* clusterLocalTrackY		= 0;
	std::vector<double>* clusterLocalTrackZ		= 0;
	std::vector<double>* clusterLocalpitch		= 0;
	std::vector<double>* clusterSensorThickness	= 0;
	std::vector<bool>*clusterSaturation		= 0;
	std::vector<bool>*clusterOverlapping		= 0;
	std::vector<bool>*clusterFarfromedge		= 0;
	std::vector<double>* clusterPath		= 0;
	std::vector<double>* clusterBdotX		= 0;
	std::vector<double>* clusterBdotY		= 0;
	std::vector<double>* clusterBdotZ		= 0;
	std::vector<double>* clusterBdotMag		= 0;
	std::vector<unsigned int>*clusterEvent		= 0;
	std::vector<unsigned int>*clusterStripEvent	= 0;

	tree->SetBranchAddress("clusterStripIdx",&clusterStripIdx);
	tree->SetBranchAddress("clusterAmplitudes",&clusterAmplitudes);
	tree->SetBranchAddress("clusterIdx",&clusterIdx);
	tree->SetBranchAddress("clusterCharge",&clusterCharge);
	tree->SetBranchAddress("clusterWidth",&clusterWidth);
	tree->SetBranchAddress("clusterLocalTrackPhi",&clusterLocalTrackPhi);
	tree->SetBranchAddress("clusterLocalTrackTheta",&clusterLocalTrackTheta);
	tree->SetBranchAddress("clusterSubdetid",&clusterSubdetid);
	tree->SetBranchAddress("clusterLayerwheel",&clusterLayerwheel);
	tree->SetBranchAddress("clusterDetid",&clusterDetid);
	tree->SetBranchAddress("clusterLocalTrackX",&clusterLocalTrackX);
	tree->SetBranchAddress("clusterLocalTrackY",&clusterLocalTrackY);
	tree->SetBranchAddress("clusterLocalTrackZ",&clusterLocalTrackZ);
	tree->SetBranchAddress("clusterLocalpitch",&clusterLocalpitch);
	tree->SetBranchAddress("clusterSensorThickness",&clusterSensorThickness);
	tree->SetBranchAddress("clusterSaturation",&clusterSaturation);
	tree->SetBranchAddress("clusterOverlapping",&clusterOverlapping);
	tree->SetBranchAddress("clusterFarfromedge",&clusterFarfromedge);
	tree->SetBranchAddress("clusterPath",&clusterPath);
	tree->SetBranchAddress("clusterBdotX",&clusterBdotX);
	tree->SetBranchAddress("clusterBdotY",&clusterBdotY);
	tree->SetBranchAddress("clusterBdotZ",&clusterBdotZ);
	tree->SetBranchAddress("clusterBdotMag",&clusterBdotMag);
	tree->SetBranchAddress("clusterEvent",&clusterEvent);
	tree->SetBranchAddress("clusterStripEvent",&clusterStripEvent);
	//------------------------------------

	//------------------------------------ //Declarations
	TCanvas* c1=new TCanvas("c1","c1");
	TCanvas* c2=new TCanvas("c2","c2");
	TCanvas* c3=new TCanvas("c3","c3");
	TCanvas* c4=new TCanvas("c4","c4",1200,800);
	TLegend* leg1=new TLegend(0.7,0.7,0.9,0.9);
	TLegend* leg2=new TLegend(0.7,0.7,0.9,0.9);
	TLegend* leg3=new TLegend(0.7,0.7,0.9,0.9);
	TLegend* leg4=new TLegend(0.7,0.7,0.9,0.9);
	// TLegend* leg5=new TLegend(0.75,0.75,0.9,0.9);

	TRandom3* r=new TRandom3();

	TH1F* hdEdx=new TH1F("hdEdx","",200,0,1000);

	double range[]={0,1.001,2.001,3.001,4.001,5.001,6.001,7.001,8.001};

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

	TH1F* hObsCharge=new TH1F("hObsCharge","",200,0,400);
	TH1F* hUnfoldCharge=new TH1F("hUnfoldCharge","",200,0,400);
	TH1F* hRecCharge=new TH1F("hRecCharge","",200,0,400);
	TH1F* hRecSelCharge=new TH1F("hRecSelCharge","",200,0,400);

	TH1F* hRelDiffObs=new TH1F("hRelDiffObs","",300,-2,6);
	TH1F* hRelDiffUnfold=new TH1F("hRelDiffUnfold","",300,-2,6);
	TH1F* hRelDiffRec=new TH1F("hRelDiffRec","",300,-2,6);
	TH1F* hRelDiffRecSel=new TH1F("hRelDiffRecSel","",300,-2,6);
	TH1F* hRelDiffWeight=new TH1F("hRelDiffWeight","",300,-2,6);

	TH1F* hRelDiff=new TH1F("hRelDiff","",300,-2,6);
	TH1F* hRelDiff1=new TH1F("hRelDiff1","",300,-2,6);
	TH1F* hRelDiff12=new TH1F("hRelDiff12","",300,-2,6);
	TH1F* hRelDiff2=new TH1F("hRelDiff2","",300,-2,6);
	TH1F* hRelDiff23=new TH1F("hRelDiff23","",300,-2,6);
	TH1F* hRelDiff3=new TH1F("hRelDiff3","",300,-2,6);

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

	TProfile* pUnderCluster=new TProfile("pSubCluster","Sub-cluster",16,0,15,-50,700);
	TProfile* pClusterSel=new TProfile("pClusterSel","Cleaned cluster",16,0,15,-50,700);
	TProfile* ligneThreshold=new TProfile("","",1,0,1);
	TProfile* ligneQStrip=new TProfile("","",1,0,1);


	Estimator Estim;

	double x0=0, x1=0, x2=0;
	double s0=0, s1=0;//, s2;

	unsigned int Idx=0;

	double ObsCW=0.;
	double ExpCW=0.;
	double UnfoldCW=0.;
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

	double QStrip=0.;
	double qstripfromtreatment=0.;

	std::vector<double> q;
	std::vector<double> QUnfold;
	std::vector<double> QUnfold_P;
	std::vector<double> QUnfold_M;
	std::vector<double> QUnfoldProb;

	std::vector<std::vector<double>> Cluster;
	std::vector<double> UnderCluster;

	//------------------------------------

	//------------------------------------ RelUncertaintyCrossTalkNoise
	// double UncertQStrip1=0.581025; //1%
	// double UncertQStrip5=0.378116; //5%
	double UncertQStrip10=0.268006; //10%
	// double UncertQStrip20=0.139197; //20%
	double UncertQS=0.;
	//------------------------------------


	for (unsigned int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);
		if(i%1000==0) cout<<"Event "<<i<<endl;

		for(unsigned int c=0;c<clusterWidth->size();c++){

			if(clusterSubdetid->at(c)==TIBid){
				if(clusterLayerwheel->at(c)<=2){
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

			if(clusterSubdetid->at(c)==TIDid){
				if(clusterLayerwheel->at(c) == 1){
					x0=x[4][0];
					x1=x[4][1];
					x2=x[4][2];
				}
				else if(clusterLayerwheel->at(c) == 2){
					x0=x[5][0];
					x1=x[5][1];
					x2=x[5][2];
				}
				else if(clusterLayerwheel->at(c) == 3){
					x0=x[6][0];
					x1=x[6][1];
					x2=x[6][2];
				}
			}

			if(clusterSubdetid->at(c)==TOBid){
				if(clusterLayerwheel->at(c)<=4){
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

			if(clusterSubdetid->at(c)==TECid){
				if(clusterLayerwheel->at(c) == 1){
					x0=x[7][0];
					x1=x[7][1];
					x2=x[7][2];
				}
				else if(clusterLayerwheel->at(c) == 2){
					x0=x[8][0];
					x1=x[8][1];
					x2=x[8][2];
				}
				else if(clusterLayerwheel->at(c) == 3){
					x0=x[9][0];
					x1=x[9][1];
					x2=x[9][2];
				}
				else if(clusterLayerwheel->at(c) == 4){
					x0=x[10][0];
					x1=x[10][1];
					x2=x[10][2];
				}
				else if(clusterLayerwheel->at(c) == 5){
					x0=x[11][0];
					x1=x[11][1];
					x2=x[11][2];
				}
				else if(clusterLayerwheel->at(c) == 6){
					x0=x[12][0];
					x1=x[12][1];
					x2=x[12][2];
				}
				else if(clusterLayerwheel->at(c) == 7){
					x0=x[13][0];
					x1=x[13][1];
					x2=x[13][2];
				}
			}

			if(clusterOverlapping->at(c)==0 && clusterFarfromedge->at(c)==1 && clusterSaturation->at(c)==0){

				ObsCW=clusterWidth->at(c);
				ExpCW=abs((tan(0.02)+cos(clusterLocalTrackPhi->at(c))*tan(clusterLocalTrackTheta->at(c)))*clusterSensorThickness->at(c)/clusterLocalpitch->at(c)); //We'll need Lorentz angle Theta_L 
				Idx=clusterIdx->at(c);

				for(unsigned int z=0;z<clusterStripIdx->size();z++){
					if(clusterStripIdx->at(z) == Idx){
						q.push_back(clusterAmplitudes->at(z)); 
					}
				}

				ObsCharge=Sum(q);

				double s1b=UncertaintyCrossTalk(x0,s0,x1,x2);

				QUnfold=Unfold(q,x0,x1,x2);

				UnfoldCharge=Sum(QUnfold);

				QUnfold_M=Unfold(q,x0+s0,x1+s1b,x2);

				QUnfold_P=Unfold(q,x0-s0,x1-s1b-2*s0,x2);

				double x0prob=r->Gaus(x0,2*s0); //cout<<"x0 : "<<x0<<endl;

				double x1prob=r->Gaus(x1,2*s1); //cout<<"x1 : "<<x1<<endl;

				double x2prob=(double)(1-2*x1-x0)/2.; //cout<<"x2 : "<<x2<<endl;

				QUnfoldProb=Unfold(q,x0prob,x1prob,x2prob); //for(unsigned int i=0;i<QUnfoldProb.size();i++) cout<<"Q : "<<QUnfoldProb[i]<<endl; getchar();

				Cluster=Cutting(Threshold(QUnfold,clusterThreshold));

				if(Cluster.size()==0) cerr<<"Warning"<<"\tEvent "<<i<<"\tCluster"<<c<<endl;
				else UnderCluster=SelectionAfterCutting(Cluster);

				for(unsigned int i=0;i<Cluster.size();i++){ //UnfoldClusterWidth
					UnfoldCW += Cluster[i].size();
				}

				double k=Cluster.size();
				if (k>1 && k<10) UnfoldCW+=k-1;
//				for(unsigned int k=2;k<10;k++){
//					if(Cluster.size()==k) UnfoldCW+=(k-1);
//				}

				QStrip=funcQStrip(UnderCluster,clusterPath,c);

				RecCW=TotalClusterWidth(UnderCluster,QStrip)[1][0];
				RecCharge=TotalClusterWidth(UnderCluster,QStrip)[1][1];

				UncertQS=UncertQStrip10;

				RecSelCW=TotalClusterWidthAfterTreatment(QUnfold,UncertQS,clusterPath,c,clusterThreshold)[0];
				RecSelCWUncert_P=TotalClusterWidthAfterTreatment(QUnfold_P,UncertQS,clusterPath,c,clusterThreshold)[0];
				RecSelCWUncert_M=TotalClusterWidthAfterTreatment(QUnfold_M,UncertQS,clusterPath,c,clusterThreshold)[0];
				RecSelCharge=TotalClusterWidthAfterTreatment(QUnfold,UncertQS,clusterPath,c,clusterThreshold)[1];

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

// new?

				if(ObsCW ==8 && ExpCW <= 2 && RecSelCW <= 2){
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
					for(unsigned int l=0;l<q.size()+2;l++) ligneThreshold->Fill(l,13.);
					pUnfoldCluster->Reset();
					for(unsigned int l=0;l<QUnfold.size();l++){
						if(QUnfold[l]<=clusterThreshold) QUnfold[l]=0.;
						pUnfoldCluster->Fill(l+1,QUnfold[l]);
					}
					for(unsigned int l=0;l<UnderCluster.size();l++){
						pUnderCluster->Fill(l+1,UnderCluster[l]);
					}

					for(unsigned int l=0;l<UnderCluster.size()+2;l++) ligneQStrip->Fill(l,qstripfromtreatment);
					for(unsigned int l=0;l<ClusterSel.size();l++){
						pClusterSel->Fill(l+1,ClusterSel[l]);
					}


					TLegend* legend=new TLegend(0.75,0.75,0.9,0.9);
					c4->cd();
					c4->Divide(2,2);
					gStyle->SetOptStat(0);
					c4->cd(1);
					pObsCluster->GetXaxis()->SetTitle("Strip");
					pObsCluster->GetYaxis()->SetTitle("Charge");
					pObsCluster->GetYaxis()->SetRangeUser(0,pUnfoldCluster->GetMaximum()+9);
					pObsCluster->SetLineColor(4);
					pObsCluster->Draw();
					pUnfoldCluster->SetLineColor(2);
					pUnfoldCluster->Draw("SAME");
					ligneThreshold->SetLineColor(1);
					ligneThreshold->SetLineStyle(7);
					ligneThreshold->Draw("SAME");
					legend->AddEntry(pObsCluster,"Observed","l");
					legend->AddEntry(pUnfoldCluster,"Unfold","l");
					legend->Draw("SAME");
//					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.pdf");
//					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.png");
					c4->cd(2);
					pUnfoldCluster->GetXaxis()->SetTitle("Strip");
					pUnfoldCluster->GetYaxis()->SetTitle("Charge");
					pUnfoldCluster->Draw();
//					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.pdf");
//					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.png");
					c4->cd(3);
					pUnderCluster->GetXaxis()->SetTitle("Strip");
					pUnderCluster->GetYaxis()->SetTitle("Charge");
					pUnderCluster->Draw();
					ligneQStrip->SetLineColor(1);
					ligneQStrip->SetLineStyle(7);
					// if(UnderCluster.size()>=4) 
					ligneQStrip->Draw("SAME");
//					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.pdf");
//					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.png");
					c4->cd(4);
					pClusterSel->GetXaxis()->SetTitle("Strip");
					pClusterSel->GetYaxis()->SetTitle("Charge");
					pClusterSel->Draw();
					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.pdf");
					c4->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/chargeCluster.png");
					getchar();
					pObsCluster->Reset();
					pUnfoldCluster->Reset();
					pUnderCluster->Reset();
					pClusterSel->Reset();
					ligneThreshold->Reset();
					ligneQStrip->Reset();
					c4->Clear();

				}

// end new?


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
				hUnfoldCharge->Fill(UnfoldCharge);
				hRecCharge->Fill(RecCharge);
				hRecSelCharge->Fill(RecSelCharge);
				hRelDiffCrossTalkWeighted->Fill(RelDiffCrossTalkWeighted);
				hdEdx->Fill(clusterCharge->at(c)/clusterPath->at(c));
				hRelDiffObs->Fill(RelDiffObs);
				hRelDiffUnfold->Fill(RelDiffUnfold);
				hRelDiffRec->Fill(RelDiffRec);
				hRelDiffRecSel->Fill(RelDiffRecSel);
				hRelDiffWeight->Fill(RelDiffWeight);

				q.clear();
				QUnfold.clear();
				QUnfold_P.clear();
				QUnfold_M.clear();
				QUnfoldProb.clear();
				Cluster.clear();
				UnderCluster.clear();
				ExpCW=0.;
				UnfoldCW=0;
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

	TFile ofile(outfilename.c_str(),"RECREATE");

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

	pObsCluster->Write();
	pUnfoldCluster->Write();


	hExpVsRec->GetXaxis()->SetTitle("Expected Cluster Width");
	hExpVsRec->GetYaxis()->SetTitle("Reconstructed Cluster Width");
	hExpVsRec->Write();
	c2->cd();
	gStyle->SetOptStat(0);
	hExpVsRec->Draw();
	hExpVsRec->SetDrawOption("COLZ");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsRec.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsRec.png");

	hExpVsObs->GetXaxis()->SetTitle("Expected Cluster Width");
	hExpVsObs->GetYaxis()->SetTitle("Observed Cluster Width");
	hExpVsObs->Write();
	c2->cd();
	gStyle->SetOptStat(0);
	hExpVsObs->Draw();
	hExpVsObs->SetDrawOption("COLZ");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsObs.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsObs.png");

	hExpVsRecSel->GetXaxis()->SetTitle("Expected Cluster Width");
	hExpVsRecSel->GetYaxis()->SetTitle("Reconstructed and Selected Cluster Width");
	hExpVsRecSel->Write();
	c2->cd();
	gStyle->SetOptStat(0);
	hExpVsRecSel->Draw();
	hExpVsRecSel->SetDrawOption("COLZ");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsRecSel.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsRecSel.png");

	hExpVsUnfold->GetXaxis()->SetTitle("Expected Cluster Width");
	hExpVsUnfold->GetYaxis()->SetTitle("Unfolded Cluster Width");
	hExpVsUnfold->Write();
	c2->cd();
	gStyle->SetOptStat(0);
	hExpVsUnfold->Draw();
	hExpVsUnfold->SetDrawOption("COLZ");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsUnfold.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsUnfold.png");

	hExpVsWeight->GetXaxis()->SetTitle("Expected Cluster Width");
	hExpVsWeight->GetYaxis()->SetTitle("Weighted Cluster Width");
	hExpVsWeight->Write();
	c2->cd();
	gStyle->SetOptStat(0);
	hExpVsWeight->Draw();
	hExpVsWeight->SetDrawOption("COLZ");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsWeight.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ExpVsWeight.png");

	hRelDiffCrossTalkWeighted->GetXaxis()->SetTitle("(Weighted - Normal)/Normal CrossTalk Cluster Width");
	c2->cd();
	gStyle->SetOptStat(0);
	hRelDiffCrossTalkWeighted->Draw();
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiffWeighted.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiffWeighted.png");


	//--------- //ClusterWidth

	c1->cd();
	c1->SetLogy();
	hExpCW->GetXaxis()->SetTitle("Cluster Width");
	hExpCW->SetLineColor(38);
	hObsCW->SetLineColor(1);
	hUnfoldCW->SetLineColor(42);
	hRecCW->SetLineColor(46);
	hRecSelCW->SetLineColor(8);
	hWeightedWithUncertaintiesCrossTalkCW->SetLineColor(40);
	hRecSelProbCW->SetLineColor(28);
	hExpCW->Draw();
	hObsCW->Draw("SAME");
	hUnfoldCW->Draw("SAME");
	hRecCW->Draw("SAME");
	hRecSelCW->Draw("SAME");
	hWeightedWithUncertaintiesCrossTalkCW->Draw("SAME");
	// hRecSelProbCW->Draw("SAME");
	leg1->AddEntry(hObsCW,"Obs","l");
	leg1->AddEntry(hExpCW,"Exp","l");
	leg1->AddEntry(hUnfoldCW,"Unfold","l");
	leg1->AddEntry(hRecCW,"Rec","l");
	leg1->AddEntry(hRecSelCW,"RecSel","l");
	leg1->AddEntry(hWeightedWithUncertaintiesCrossTalkCW,"WeightedCrossTalk","l");
	// leg1->AddEntry(hRecSelProbCW,"RecSelProb","l");
	leg1->Draw("SAME");
	c1->Write();
	gStyle->SetOptStat(0);
	c1->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ClusterWidth.pdf");
	c1->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/ClusterWidth.png");

	//---------

	//--------- //RelDiffCWTopo

	c2->cd();
	hRelDiff2->GetXaxis()->SetTitle("(Rec-Exp)/Exp Cluster Width");
	hRelDiff1->SetLineColor(1);
	hRelDiff12->SetLineColor(38);
	hRelDiff2->SetLineColor(42);
	hRelDiff23->SetLineColor(46);
	hRelDiff3->SetLineColor(8);
	hRelDiff1Weighted->SetLineColor(1);
	hRelDiff12Weighted->SetLineColor(38);
	hRelDiff2Weighted->SetLineColor(42);
	hRelDiff23Weighted->SetLineColor(46);
	hRelDiff3Weighted->SetLineColor(8);
	// hRelDiff1Weighted->SetLineStyle(2);
	// hRelDiff12Weighted->SetLineStyle(2);
	// hRelDiff2Weighted->SetLineStyle(2);
	// hRelDiff23Weighted->SetLineStyle(2);
	// hRelDiff3Weighted->SetLineStyle(2);
	hRelDiff2->Draw();
	hRelDiff1->Draw("SAME");
	hRelDiff12->Draw("SAME");
	hRelDiff23->Draw("SAME");
	hRelDiff3->Draw("SAME");
	leg2->AddEntry(hRelDiff1,"CW=1","l");
	leg2->AddEntry(hRelDiff12,"1<CW<2","l");
	leg2->AddEntry(hRelDiff2,"CW=2","l");
	leg2->AddEntry(hRelDiff23,"2<CW<3","l");
	leg2->AddEntry(hRelDiff3,"CW>=3","l");
	leg2->Draw("SAME");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiff.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiff.png");
	c2->Write();
	c2->cd();
	hRelDiff1Weighted->Draw();
	hRelDiff12Weighted->Draw("SAME");
	hRelDiff2Weighted->Draw("SAME");
	hRelDiff23Weighted->Draw("SAME");
	hRelDiff3Weighted->Draw("SAME");
	leg2->Draw("SAME");
	c2->Write();
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiff2.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiff2.png");



	//---------

	//--------- //RelDiffCrossTalk

	c2->cd();
	hRelDiffCrossTalk_M->GetXaxis()->SetTitle("(RecWithDiffCrossTalk-RecSel)/RecSel Cluster Width");
	hRelDiffCrossTalk_P->SetLineColor(1);
	hRelDiffCrossTalk_M->SetLineColor(42);
	// hRelDiff->SetLineColor(46);
	hRelDiffCrossTalkWeighted->SetLineColor(38);
	// hRelDiffCrossTalk_P->Scale(1./hRelDiffCrossTalk_P->Integral());
	// hRelDiffCrossTalk_M->Scale(1./hRelDiffCrossTalk_M->Integral());
	hRelDiffCrossTalk_M->Draw();
	hRelDiffCrossTalk_P->Draw("SAME");
	// hRelDiff->Draw("SAME");
	hRelDiffCrossTalkWeighted->Draw("SAME");

	// leg3->AddEntry(hRelDiff,"Rec","l");
	leg3->AddEntry(hRelDiffCrossTalk_P,"CrossTalk_P","l");
	leg3->AddEntry(hRelDiffCrossTalk_M,"CrossTalk_M","l");
	leg3->AddEntry(hRelDiffCrossTalkWeighted,"CrossTalkWeighted","l");
	leg3->Draw("SAME");
	c2->Write();
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiffCrossTalk.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiffCrossTalk.png");

	//---------

	//--------- //Charge
	c3->cd();
	hObsCharge->SetLineColor(1);
	hUnfoldCharge->SetLineColor(38);
	hRecCharge->SetLineColor(46);
	hRecSelCharge->SetLineColor(8);
	hObsCharge->Draw();
	hUnfoldCharge->Draw("SAME");
	hRecCharge->Draw("SAME");
	hRecSelCharge->Draw("SAME");
	leg4->AddEntry(hObsCharge,"Obs","l");
	leg4->AddEntry(hUnfoldCharge,"Unfold","l");
	leg4->AddEntry(hRecCharge,"Rec","l");
	leg4->AddEntry(hRecSelCharge,"RecSel","l");
	leg4->Draw("SAME");
	c3->Write();
	c3->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/Charge.pdf");
	c3->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/Charge.png");
	//---------

	//--------- //RelDiff

	TLegend* legRelDiff=new TLegend(0.7,0.7,0.9,0.9);
	c2->cd();  
	hRelDiffObs->GetXaxis()->SetTitle("(CW - ExpCW)/ExpCW");
	hRelDiffObs->GetYaxis()->SetRangeUser(0,hRelDiffRec->GetMaximum()+9);
	hRelDiffObs->SetLineColor(1);
	hRelDiffUnfold->SetLineColor(42);
	hRelDiffRec->SetLineColor(46);
	hRelDiffRecSel->SetLineColor(8);
	hRelDiffWeight->SetLineColor(40);
	hRelDiffObs->Draw();
	hRelDiffUnfold->Draw("SAME");
	hRelDiffRec->Draw("SAME");
	hRelDiffRecSel->Draw("SAME");
	hRelDiffWeight->Draw("SAME");
	legRelDiff->AddEntry(hRelDiffObs,"Obs","l");
	legRelDiff->AddEntry(hRelDiffUnfold,"Unfold","l");
	legRelDiff->AddEntry(hRelDiffRec,"Rec","l");
	legRelDiff->AddEntry(hRelDiffRecSel,"RecSel","l");
	legRelDiff->AddEntry(hRelDiffWeight,"Weight","l");
	legRelDiff->Draw("SAME");
	c2->Write();
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiff2.pdf");
	c2->SaveAs("/afs/hephy.at/user/l/llechner/www/ClusterSplitting/MCbased_v1/RelDiff2.png");



	//---------


	ofile.Write();
	ofile.Close();

	return 0;
}
