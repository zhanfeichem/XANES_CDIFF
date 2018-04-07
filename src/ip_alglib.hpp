#ifndef IP_LIB
#define IP_LIB


#include <iostream>
#include <fstream>
#include<sstream>
#include<vector>
#include<math.h>
 
#include <stdafx.h>
#include <interpolation.h>
 
#include<stdio.h>
#include<stdlib.h>

using namespace std;
using namespace alglib;


int ip(vector<double> x1,vector<double> y1,vector<double> x2,vector<double> y2,vector<double>& xr,vector<double>& yr1,vector<double>& yr2,int np);

int linspace(double a,double b,int n,vector<double> & vec_out);
int getint(vector<double> x1,vector<double> x2,double& a,double& b);


int ip(vector<double> x1,vector<double> y1,vector<double> x2,vector<double> y2,vector<double>& xr,vector<double>& yr1,vector<double>& yr2,int np=1000)
{
     
     
     
    double* ptr_yn1;
    double* ptr_yn2;
    ofstream out1("out.txt");
    real_1d_array xx1;
    real_1d_array yy1;
    real_1d_array xx2;
    real_1d_array yy2;
    real_1d_array xxn;
    real_1d_array yy1n;
    real_1d_array yy2n;
    real_1d_array dd1n;
    real_1d_array dd2n;
    vector<double> xn_vec; 
    double anew,bnew = 0;
    getint(x1,x2,anew,bnew);
    linspace(anew,bnew,np,xn_vec);
    xxn.setcontent(xn_vec.size(),&xn_vec[0]);
    xx1.setcontent(x1.size(),&x1[0]);
    yy1.setcontent(y1.size(),&y1[0]);
    xx2.setcontent(x2.size(),&x2[0]);
    yy2.setcontent(y2.size(),&y2[0]);
    spline1dconvdiffcubic(xx1,yy1,xxn,yy1n,dd1n);
    spline1dconvdiffcubic(xx2,yy2,xxn,yy2n,dd2n);
    ptr_yn1=yy1n.getcontent();
    ptr_yn2=yy2n.getcontent();
    vector<double> yn1_vec(ptr_yn1,ptr_yn1+yy1n.length()); 
    vector<double> yn2_vec(ptr_yn2,ptr_yn2+yy2n.length()); 

     
    xr.clear();
    yr1.clear();
    yr2.clear();
    xr=xn_vec;
    yr1=yn1_vec;
    yr2=yn2_vec;

} 
int getint(vector<double> x1,vector<double> x2,double& a,double& b){
    if(*x1.begin()>=*x2.begin()) a=*x1.begin();else a=*x2.begin();
    if(*(x1.end()-1)<=*(x2.end()-1)) b=*(x1.end()-1); else b=*(x2.end()-1);
    if(a>=b) {cerr<<"getint error"<<endl;exit(-1);}
    return 0;
}
int linspace(double a,double b,int n,vector<double> & vec_out){
    vec_out.clear();
    vec_out.push_back(a);
    double h=(b-a)/(n-1);
    for(int i=1;i<=n-2;++i){
        double cc=a+i*h;
        vec_out.push_back(cc);
    } 
    vec_out.push_back(b);
    return 0;
} 


#endif  
