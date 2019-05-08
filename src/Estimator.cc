#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <algorithm>

#include "../interface/Estimator.h"
//#include "ClusterShape/interface/Estimator.h"

using namespace std;

Estimator::Estimator(){
	Vect_.push_back(1);
}

Estimator::Estimator(vector<double> Vect){
	Vect_ = Vect;
}

Estimator::~Estimator(){
	Vect_.clear();
}

void Estimator::setVect(vector<double> Vect){
	Vect_ = Vect;
}

vector<double> Estimator::getVect(){
	return Vect_;
}

void Estimator::dispVect(){
	for(unsigned int i=0;i<Vect_.size();i++){cout<<"Element "<<i<<" du vecteur : "<<Vect_[i]<<endl;}
}

double Estimator::mean(){
	double mean = .0;
	unsigned int size = Vect_.size();
	for(unsigned int i=0;i<size;i++){
		mean+=Vect_[i];}
	return mean/size;
}

double Estimator::mean(vector<double> Vect){
	double mean = .0;
	unsigned int size = Vect.size();
	for(unsigned int i=0;i<size;i++){
		mean+=Vect[i];}
	return mean/size;
}

double Estimator::median(){
	sort(Vect_.begin(),Vect_.end());
	unsigned int size = Vect_.size();
	double med = Vect_[size/2.];
	return med;}

double Estimator::trunc40(){
	sort(Vect_.begin(),Vect_.end());
	unsigned int size = Vect_.size();
	double size60 = size*0.6;
	vector<double> Vect60;
	for(unsigned int i=0; i<size60; i++){Vect60.push_back(Vect_[i]);}
	return mean(Vect60);}

double Estimator::trunc40weight(vector<double> vpath){
    sort(Vect_.begin(),Vect_.end());
    unsigned int size = Vect_.size();
    double size60 = size*0.6;
    vector<double> Vect60;
    vector<double> path60;
    for(unsigned int i=0; i<size60; i++){
		Vect60.push_back(Vect_[i]);
        path60.push_back(vpath[i]);}
    double res=0;
    double w=0;
    double wtot=0;
    for(unsigned int i=0; i<size60; i++){
        w=pow(path60[i], 0.5);
        res+=w*Vect60[i];
        wtot+=w;}
    res=res/wtot;
    return res;}


double Estimator::harmonic2(){
	double mean = .0;
	unsigned int size = Vect_.size();
	for(unsigned int i=0;i<size;i++){
		mean+=pow(Vect_[i],-2);
	}
	return pow(mean/size,-0.5);}

double Estimator::stdDev(){
	unsigned int size = Vect_.size();
	double res = .0;
	for(unsigned int i=0; i<size; i++){
	res+=pow(Vect_[i]-mean(),2);}
	res=pow(res/(double)(size-1),0.5);
	return res;}

double Estimator::weight(vector<double> vpath){
    unsigned int size = Vect_.size();
    double res = 0;
    double w = 0;
    double wtot = 0;
    for(unsigned int i=0; i<size; i++){
    	w = pow(vpath[i]/0.32,0.266);
    	res += w*Vect_[i];
    	wtot += w;
    }
    res=res/wtot;
    return res;}

double Estimator::trunchl(){
    sort(Vect_.begin(),Vect_.end());
    unsigned int size = Vect_.size();
    double size1 = size*0.15;
    double size2 = size*(1.0-0.35);
    vector<double> Vect1;
    for(int i=size1; i<size2; i++){Vect1.push_back(Vect_[i]);}
    return mean(Vect1);}

double Estimator::meanWithoutFL(){
	unsigned int size = Vect_.size();
	unsigned int size2 = size-1;
	vector<double> Vect2;
	for(unsigned int i=1;i<size2;i++){Vect2.push_back(Vect_[i]);}
	return mean(Vect2);
}

int Estimator::posMax(){
	unsigned int size = Vect_.size();
	double max=Vect_[0];
	int index=0;
	for(unsigned int i=0;i<size;i++){
		if(Vect_[i]>max) max=Vect_[i], index=i;
	}
	return index;
}

int Estimator::posMin(){
	unsigned int size = Vect_.size();
	double min=Vect_[0];
	unsigned int index=0;
	for(unsigned int i=0;i<size;i++){
		if(Vect_[i]<min) min=Vect_[i], index=i;
	}
	return index;
}

int Estimator::posMin0(){
	unsigned int size = Vect_.size();
	double min=Vect_[0];
	unsigned int index=0;
	for(unsigned int i=0;i<size;i++){
		if(Vect_[i]<min && Vect_[i]!=0) min=Vect_[i], index=i;
	}
	return index;
}
