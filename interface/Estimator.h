#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <algorithm>

using namespace std;

class Estimator{

private :

 vector<double> Vect_;

public :

 Estimator();
 Estimator(vector<double>);
 ~Estimator();
 void setVect(vector<double>);
 vector<double> getVect();
 void dispVect();
 double mean();
 double mean(vector<double>);
 double median();
 double trunc40();
 double harmonic2();
 double stdDev();
 double weight(vector<double>);
 double trunchl();
 double trunc40weight(vector<double>);
 double meanWithoutFL();
 double dEdx_PredWidth(double);
 int posMax();
 int posMin();
 int posMin0();
};
