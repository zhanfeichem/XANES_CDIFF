#ifndef IP_LIB
#define IP_LIB

#include <iostream>
#include <vector>
#include <fstream>

#define DIS
using namespace std;
int ip(vector<double> x1,vector<double> y1,vector<double> x2,vector<double> y2,vector<double>& xr,vector<double>& yr1,vector<double>& yr2,int np);
int ip3(vector<double> x1,vector<double> y1,vector<double> x2,vector<double> y2,vector<double> x3,vector<double> y3,vector<double>& xr,vector<double>& yr1,vector<double>& yr2,vector<double>& yr3,int np);
int getint(vector<double> x1,vector<double> x2,double& a,double& b);
int linspace(double a,double b,int n,vector<double> & vec_out);

int interp1d(vector<double> x,vector<double> y,vector<double> xn,vector<double>& yn);
void spline(vector<double> &x, vector<double> &y, const double yp1, const double ypn,vector<double> &y2);
void splint(vector<double> &xa, vector<double> &ya, vector<double> &y2a, const double x, double &y);


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


 
 
 
 
 
 

int ip(vector<double> x1,vector<double> y1,vector<double> x2,vector<double> y2,vector<double>& xr,vector<double>& yr1,vector<double>& yr2,int np=1000){



xr.clear();
yr1.clear();
yr2.clear();
double anew=0;double bnew=0;
getint(x1,x2,anew,bnew);
linspace(anew,bnew,np,xr);
interp1d(x1,y1,xr,yr1);
interp1d(x2,y2,xr,yr2);
return 0;
} 

int ip3(vector<double> x1,vector<double> y1,vector<double> x2,vector<double> y2,vector<double> x3,vector<double> y3,vector<double>& xr,vector<double>& yr1,vector<double>& yr2,vector<double>& yr3,int np=1000){
 
xr.clear();
yr1.clear();
yr2.clear();
yr3.clear();
vector<double> xr_int;
double anew=0;double bnew=0;
getint(x1,x2,anew,bnew);
linspace(anew,bnew,np,xr_int);
getint(xr_int,x3,anew,bnew); 
linspace(anew,bnew,np,xr);
interp1d(x1,y1,xr,yr1);
interp1d(x2,y2,xr,yr2);
interp1d(x3,y3,xr,yr3);
return 0;
} 

int interp1d(vector<double> x,vector<double> y,vector<double> xn,vector<double>& yn){
     
    #ifdef DIS
    if(x.size()!=y.size()) {cerr<<"interp1d error x y do not match"<<endl;exit(-1);}
    #endif  
    yn.clear();
    vector<double> der;
    spline(x,y,0,0,der);

    double y_val=0;
    for(unsigned i=0;i<xn.size();++i){
            splint(x,y,der,xn[i],y_val);
            yn.push_back(y_val);

             
    }

return 0;} 

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



void spline(vector<double> &x, vector<double> &y, const double yp1, const double ypn,vector<double> &y2){
	int i,k;
	double p,qn,sig,un;

	int n=y.size();
	vector<double> u(n,0); 
	y2.clear();
	y2=u; 
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
} 


void splint(vector<double> &xa, vector<double> &ya, vector<double> &y2a, const double x, double &y){
     
    int k;
	double h,b,a;
	int n=xa.size();
	 
    if(x>=xa[n-1]){y=ya[n-1];return ;}
    if(x<=xa[0]){y=ya[0];return;}
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	 
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
} 




#endif  
