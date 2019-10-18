#include <iostream>
#include <vector>
#include <list>
#include <limits>
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

#include "ClusterShape/interface/Estimator.h"

using namespace std;

vector<double> Unfold(const vector<double>& q, double x0, double x1, double x2){
    // Create a SVector representing the charges
    constexpr int Nsize=12;
    ROOT::Math::SVector<double, Nsize> v;
    unsigned int max=q.size();
    if(q.size()>Nsize) max=Nsize;
    for(unsigned int i=0;i<Nsize;i++){
        v(i)=0;
    }
    for(unsigned int i=0;i<max;i++){
        v(i)=q[i];
    }

    // Create a matrix for the cross-talk
    ROOT::Math::SMatrix<double,Nsize> m;

    for(unsigned int i=0;i<Nsize;i++){
        m(i,i)=x0;
        if(i<(Nsize-1)){
            m(i,i+1)=x1;
            m(i+1,i)=x1;
        }

        if(i<(Nsize-2)){
            m(i,i+2)=x2;
            m(i+2,i)=x2;
        }
    }

    // Invert the matrix
    m.Invert();

    // Apply the inverted matrix on the cluster amplitude
    ROOT::Math::SVector<double,Nsize> vOut=v*m;
    vector<double> out(q.size());
    for(unsigned int i=0;i<q.size();i++){
        out[i]=vOut(i);
    }

    return out;
}

double QperStrip(double path, double Q1, double Q2, int n){
    double medQdist=300; // change that to be flexible?
    return (medQdist*path-(Q1+Q2))/n;
}

double min(double a, double b){
    return (a < b) ? a : b;
}

double Mean(vector<double> Vect){
    double mean=.0;
    unsigned int size=Vect.size();
    for(unsigned int i=0;i<size;i++){
        mean+=Vect[i];
    }
    return mean/size;
}

int PosMax(vector<double> Vect_){
    unsigned int size=Vect_.size();
    double max=Vect_[0];
    unsigned int index=0;
    for(unsigned int i=0;i<size;i++){
        if(Vect_[i]>max) max=Vect_[i], index=i;
    }
    return index;
}

double MeanHarmonic(vector<double> Vect){
    double mean=0.;
    for(unsigned int i=0;i<Vect.size();i++){
        mean+=pow(Vect[i],2);
    }
    return sqrt(mean);
}

double Sum(vector<double> Vect){
    double res=0;
    for(unsigned int i=0;i<Vect.size();i++){
        res+=Vect[i];
    }
    return res;
}

double UncertaintyCrossTalk(double x0, double s0, double x1, double x2){
    return (double)(1-x0-s0)/2.-x1-x2;
}

vector<double> Threshold(vector<double> Q, double clusterThreshold){
    // return charge vector where charges below threshold are set to 0 (for later cutting)
    std::vector<double> QThreshold;
    for(unsigned int i=0;i<Q.size();i++){
        if(Q[i]>clusterThreshold){
            QThreshold.push_back(Q[i]);
        }
        else{
            QThreshold.push_back(0);
        }
    }
    return QThreshold;
}

double funcQStrip(vector<double> Cluster, vector<double>* clusterPath, int c){
    Estimator Estim;
    double qstrip=0.;

    if(Cluster.size()>2){
        // QperStrip only works for >2 as it uses edge hits
        Estim.setVect(Cluster);
        qstrip=QperStrip( clusterPath->at(c), Cluster[0], Cluster[Cluster.size()-1], Cluster.size()-2 );
        double meanWithoutFL=Estim.meanWithoutFL();

        if(!(qstrip>0 && qstrip<2*meanWithoutFL)){
            qstrip=Estim.meanWithoutFL();
        }
    }

    return qstrip;
}

vector<vector<double>> Cutting(vector<double> Q){
    // remove edge entries which are 0
    // create new vector if there is a gap

    std::vector<std::vector<double>> Cluster;
    std::vector<double> UnderCluster;

    for(unsigned int k=0;k<Q.size();k++){
        if(Q[k]!=0) UnderCluster.push_back(Q[k]);
        if(Q[k]==0 && UnderCluster.size()!=0){
            Cluster.push_back(UnderCluster);
            UnderCluster.clear();
        }
        if(Q[k] != 0 && k==Q.size()-1) Cluster.push_back(UnderCluster);
    }

    UnderCluster.clear();
    return Cluster;
}

vector<double> SelectionAfterCutting(vector<vector<double>> Cluster){
    // Select cluster with max charge
    std::vector<double> Interm;
    Interm.reserve(Cluster.size());

    for(unsigned int k=0;k<Cluster.size();k++){
        Interm.push_back(Sum(Cluster[k]));
    }

    return Cluster[PosMax(Interm)];
}

vector<double> ClusterizationWithNoise(vector<double> UnderCluster, double QStrip, double RelDiffQStrip){
    // Set noise bins to 0, cut clusters and select the one with largest charge
    if(UnderCluster.size()<4){
        return UnderCluster;
    }
    else{
        double thresh=QStrip*(1-RelDiffQStrip);

        for(unsigned int i=1;i<UnderCluster.size()-1;i++){
            if(UnderCluster[i]<thresh) {
                UnderCluster[i]=0;
            }
        }

        return SelectionAfterCutting(Cutting(UnderCluster));
    }
}

vector<vector<double>> TotalClusterWidth(vector<double> Cluster, double QStrip){

    std::vector<std::vector<double>> res;
    res.reserve(2);
    std::vector<double> ClusterWidth;
    ClusterWidth.reserve(Cluster.size());

    if(Cluster.size()==1){
        ClusterWidth.push_back(1);
    }
    else if(Cluster.size()==2){
        ClusterWidth.push_back(1);
        ClusterWidth.push_back(1);
    }
    else if(Cluster.size()>2){
        // add partial charges for edge strips, 1 for intermed strips

        ClusterWidth.push_back(min(Cluster[0]/QStrip,1));
        for(int j=1;j<Cluster.size()-1;j++){
            ClusterWidth.push_back(1);
        }
        ClusterWidth.push_back(min(Cluster[Cluster.size()-1]/QStrip,1));

    }

    res.push_back(ClusterWidth);
    res.push_back({Sum(ClusterWidth),Sum(Cluster)});

    return res;
}

vector<double> TotalClusterWidthAfterTreatment(vector<double> QUnfold, double UncertQStrip, vector<double>* clusterPath, int c, double clusterThreshold){

    Estimator Estim;
    std::vector<std::vector<double>> TotalClusterWidthVec;
    std::vector<std::vector<double>> Cluster;
    std::vector<double> UnderCluster;
    std::vector<double> res;
    res.reserve(3);
    double QStrip=0.;
    double RecSelCW=0.;
    double RecSelTotCharge=0.;
    double qcut=0;

    Cluster=Cutting(Threshold(QUnfold, clusterThreshold));

    if(Cluster.size()!=0){
        UnderCluster=SelectionAfterCutting(Cluster);
    }

    QStrip=funcQStrip(UnderCluster,clusterPath,c);
    qcut=QStrip*(1-UncertQStrip);

    UnderCluster=ClusterizationWithNoise(UnderCluster,QStrip,UncertQStrip);
    QStrip=funcQStrip(UnderCluster,clusterPath,c);
    TotalClusterWidthVec=TotalClusterWidth(UnderCluster,QStrip);

    RecSelCW=TotalClusterWidthVec[1][0];
    RecSelTotCharge=TotalClusterWidthVec[1][1];

    Cluster.clear();
    UnderCluster.clear();

    res.push_back(RecSelCW);
    res.push_back(RecSelTotCharge);
    res.push_back(qcut);

    return res;
}


void list_dir(const char *path, vector<string> &stringvec){
    struct dirent *entry;
    DIR *dir=opendir(path);

    if(dir==NULL){
      return;
    }

    unsigned int i=0;
    while((entry=readdir(dir)) != NULL){
        stringvec[i]=entry->d_name;
        cout << stringvec[i];
        i++;
    }

    closedir(dir);

}
