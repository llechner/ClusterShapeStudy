#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <dirent.h>
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <sys/types.h>

int PosMax(vector<double> Vect_);
double Sum(vector<double> Vect);
double Mean(vector<double> Vect);
double MeanHarmonic(vector<double> Vect);

double QperStrip(double path, double Q1, double Q2, int n);
double funcQStrip(vector<double> Cluster, vector<double>* clusterPath, int c);
double UncertaintyCrossTalk(double x0, double s0, double x1, double x2);

vector<double> Unfold(const vector<double>& q, double x0, double x1, double x2);
vector<double> Threshold(vector<double> Q, double clusterThreshold);
vector<vector<double>> Cutting(vector<double> Q);
vector<double> SelectionAfterCutting(vector<vector<double>> Cluster);

vector<double> ClusterizationWithNoise(vector<double> UnderCluster, double QStrip, double RelDiffQStrip);
vector<vector<double>> TotalClusterWidth(vector<double> Cluster, double QStrip);
vector<double> TotalClusterWidthAfterTreatment(vector<double> QUnfold, double UncertQStrip, vector<double>* clusterPath, int c, double clusterThreshold);

void list_dir(const char *path, vector<string> &stringvec);
