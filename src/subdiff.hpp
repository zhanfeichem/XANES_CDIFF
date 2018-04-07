#ifndef SUBDIFF
#define SUBDIFF


#include <iostream>
#include <fstream>
#include<sstream>
#include<vector>

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

 
 
#include "ipmy.hpp"
#include "read_dat.hpp"
#include "nlopt.h"
#include "glo.hpp"

namespace feffdiff{ 

using namespace std;

struct opt_trans{
    vector<double> x_xmu;
    vector<double> y_xmu;
    vector<double> x_prm_the;
    vector<double> y_prm_the;
    vector<double> x_diff_exp;
    vector<double> y_diff_exp;
    double es_prm;
    double ra_prm;
}; 

double trapz(vector<double> y,double h);
double subopt(int n_fro,double* x_fro,char xmupos_char[]);
double optfun(unsigned n, const double *x, double *grad, void * f_data);
double optfun_save(int n_fro,double* x_fro,int n,const double *x,void * f_data);
 
 
 
 
double subopt(int n_fro,double* x_fro,char xmupos_char[]){ 
     


    vector<vector<double> > dat_xmu; 
    vector<vector<double> > dat_prm_the; 
    vector<vector<double> > dat_diff_exp; 
    vector<double> x_xmu,y_xmu,x_prm_the,y_prm_the,x_diff_exp,y_diff_exp;


     
    read_dat(xmupos_char,6,dat_xmu);

    read_dat(myinput.diff_exp_char,2,dat_diff_exp);
    read_dat(myinput.prm_the_char,2,dat_prm_the);




    x_xmu=dat_xmu[0];y_xmu=dat_xmu[3];
     
    ofstream exfeff("1.feff");
    for(unsigned i=0;i<=x_xmu.size()-1;++i)exfeff<<x_xmu[i]<<"  "<<y_xmu[i]<<endl;
     
    x_diff_exp=dat_diff_exp[0];y_diff_exp=dat_diff_exp[1];
    x_prm_the=dat_prm_the[0];y_prm_the=dat_prm_the[1];
     
 
 
 

 
 
 
 



 
opt_trans myt; 
myt.x_xmu=x_xmu;
myt.y_xmu=y_xmu;
myt.x_prm_the=x_prm_the;
myt.y_prm_the=y_prm_the;
myt.x_diff_exp=x_diff_exp;
myt.y_diff_exp=y_diff_exp;
myt.es_prm=myinput.es_prm; 
myt.ra_prm=myinput.ra_prm; 
 


    ofstream optsave("1save.txt");
    for(unsigned i=0;i<=x_diff_exp.size()-1;++i)
         
       optsave<<x_diff_exp[i]<<"  "<<y_diff_exp[i]<<endl;
 
int npar=2;
 
 
double lb[2] = {myinput.lo_es_diff,myinput.lo_ex_diff};
double up[2] = {myinput.up_es_diff,myinput.up_ex_diff};
double x0[2] = {myinput.ini_es_diff,myinput.ini_ex_diff};   

nlopt_opt opt;
 
if(myinput.nstr_ALG=="COBYLA") opt = nlopt_create(NLOPT_LN_COBYLA,npar);
if(myinput.nstr_ALG=="BOBYQA") opt = nlopt_create(NLOPT_LN_BOBYQA,npar);
if(myinput.nstr_ALG=="PRAXIS") opt = nlopt_create(NLOPT_LN_PRAXIS,npar);
if(myinput.nstr_ALG=="NELDERMEAD") opt = nlopt_create(NLOPT_LN_NELDERMEAD,npar);
if(myinput.nstr_ALG=="SBPLX") opt = nlopt_create(NLOPT_LN_SBPLX,npar);
 
if(myinput.nstr_ALG=="ISRES") opt = nlopt_create(NLOPT_GN_ISRES,npar);
if(myinput.nstr_ALG=="ESCH") opt = nlopt_create(NLOPT_GN_ESCH,npar);
if(myinput.nstr_ALG=="DIRECT_L") opt = nlopt_create(NLOPT_GN_DIRECT_L,npar);
if(myinput.nstr_ALG=="CRS2_LM") opt = nlopt_create(NLOPT_GN_CRS2_LM,npar);
 
nlopt_set_lower_bounds(opt,lb);
nlopt_set_upper_bounds(opt,up);
nlopt_set_min_objective(opt,optfun,&myt);
 
nlopt_set_maxeval(opt,myinput.nstr_num); 
double minf; /* the minimum objective value, upon return */
 




if (nlopt_optimize(opt, x0, &minf) < 0) {
    printf("nlopt failed!\n");
}
else {
     
}
optfun_save(n_fro,x_fro,2,x0,&myt); 
 
nlopt_destroy(opt);
return minf; 
} 











double optfun(unsigned n, const double *x, double *grad, void * f_data){
     
    double es_diff = x[0]; 
    double ex = x[1]; 
    opt_trans* tt=(opt_trans *) f_data;
     
    vector<double> x_xmu = tt->x_xmu;
    vector<double> y_xmu = tt->y_xmu;
    vector<double> x_prm_the = tt->x_prm_the;
    vector<double> y_prm_the = tt->y_prm_the;
    vector<double> x_diff_exp = tt->x_diff_exp;
    vector<double> y_diff_exp = tt->y_diff_exp;
    double es_prm = tt->es_prm;
    double ra_prm = tt->ra_prm;




     
    for(vector<double>::iterator i=x_xmu.begin();i!=x_xmu.end();++i)
        *i=*i+es_prm+es_diff; 
    for(vector<double>::iterator i=y_xmu.begin();i!=y_xmu.end();++i)
        *i=(*i)*ra_prm; 
     


    vector<double> xn,yn_xmu,yn_prm_the,yn_diff_exp;

    ip(x_xmu,y_xmu,x_prm_the,y_prm_the,xn,yn_xmu,yn_prm_the,myinput.ip_num);
     

    ip(xn,yn_xmu,x_diff_exp,y_diff_exp,xn,yn_xmu,yn_diff_exp,myinput.ip_num); 
    vector<double> diff_expthe; 
    diff_expthe.clear();



    for(unsigned i=0;i<=xn.size()-1;++i){
        double dd = 0;
        dd=ex*(yn_xmu[i]-yn_prm_the[i])-yn_diff_exp[i]; 

        dd=fabs(dd);
        diff_expthe.push_back(dd);
    } 
    double result;
     
    vector<double> xout_o3; 
    vector<double> yexpout_o3,ydiffout_o3; 
    vector<double> yintout_o3,yprmout_o3; 
    if(myinput.OBJ==3)
    { 
      
        result=0;
        for(int i=0;i<x_diff_exp.size();i++)
        {
            if(x_diff_exp[i]>myinput.OBJ_aout&&x_diff_exp[i]<myinput.OBJ_bout)
                {
                xout_o3.push_back(x_diff_exp[i]);
                yexpout_o3.push_back(y_diff_exp[i]);
                } 
        } 
        interp1d(x_prm_the,y_prm_the,xout_o3,yprmout_o3);
        interp1d(x_xmu,y_xmu,xout_o3,yintout_o3);
         
        int nfit=0;
        for(unsigned i=0;i<=xout_o3.size()-1;++i)
        {
            double yi=0;
            yi=ex*(yintout_o3[i]-yprmout_o3[i]);
            ydiffout_o3.push_back(yi);
        } 
        for(unsigned i=0;i<=xout_o3.size()-1;++i)
        {
            if(xout_o3[i]>myinput.OBJ_afit&&xout_o3[i]<myinput.OBJ_bfit)
            {
            double dd = 0;
             
            dd=ydiffout_o3[i]-yexpout_o3[i];
            dd=dd*dd;
            result=result+dd;
            nfit=nfit+1;
            }
        } 
        result=result*(double)myinput.npar/(double)nfit;
        double OBJ_eps=1;
        double eps_2=OBJ_eps*OBJ_eps;
        result=result/eps_2;
    } 

    if(myinput.OBJ==2){
        result=0;
        for(unsigned i=0;i<=xn.size()-1;++i){
            double dd = 0;
            dd=ex*(yn_xmu[i]-yn_prm_the[i])-yn_diff_exp[i]; 
            dd=dd*dd;
            result=result+dd;} 
         result=result*(double)myinput.npar/(double)myinput.ip_num;
         double eps_2=myinput.OBJ_eps*myinput.OBJ_eps;
         result=result/eps_2;} 
    if(myinput.OBJ==1){
    result=trapz(diff_expthe,xn[1]-xn[0]);
    double xlen=*(xn.end()-1)-*xn.begin();
    result=result/xlen;
    } 

    if(myinput.OBJ==11)
        {
        result=0;
        double up=0;
        double down=0;
        for(unsigned i=0;i<=xn.size()-1;++i)
            {
            double dd = 0;
            dd=ex*(yn_xmu[i]-yn_prm_the[i])-yn_diff_exp[i]; 
            dd=dd*dd;
            up=up+dd;
            down=down+(yn_diff_exp[i]*yn_diff_exp[i]);
            } 
         result=up/down;
        } 



    return result;


}


double optfun_save(int n_fro,double* x_fro,int n,const double *x,void * f_data){
     
    double es_diff = x[0]; 
    double ex = x[1]; 
    opt_trans* tt=(opt_trans *) f_data;
     
    vector<double> x_xmu = tt->x_xmu;
    vector<double> y_xmu = tt->y_xmu;
    vector<double> x_prm_the = tt->x_prm_the;
    vector<double> y_prm_the = tt->y_prm_the;
    vector<double> x_diff_exp = tt->x_diff_exp;
    vector<double> y_diff_exp = tt->y_diff_exp;
     
    ofstream exfeff2("2.feff");
    for(unsigned i=0;i<=x_xmu.size()-1;++i) exfeff2<<x_xmu[i]<<"  "<<y_xmu[i]<<endl;
     
    double es_prm = tt->es_prm;
    double ra_prm = tt->ra_prm;
     
    for(vector<double>::iterator i=x_xmu.begin();i!=x_xmu.end();++i)
        *i=*i+es_prm+es_diff; 
    for(vector<double>::iterator i=y_xmu.begin();i!=y_xmu.end();++i)
        *i=(*i)*ra_prm; 
     
    vector<double> xn,yn_xmu,yn_prm_the,yn_diff_exp;
    ip(x_xmu,y_xmu,x_prm_the,y_prm_the,xn,yn_xmu,yn_prm_the,myinput.ip_num);
    ip(xn,yn_xmu,x_diff_exp,y_diff_exp,xn,yn_xmu,yn_diff_exp,myinput.ip_num); 
     
    ofstream exfeff3("3.feff");
    for(unsigned i=0;i<=x_xmu.size()-1;++i) exfeff3<<x_xmu[i]<<"  "<<y_xmu[i]<<endl;
     
    vector<double> diff_expthe; 
    vector<double> y_diff_the;
    diff_expthe.clear();y_diff_the.clear();

    for(unsigned i=0;i<=xn.size()-1;++i){
        double dd = 0;
        dd=ex*(yn_xmu[i]-yn_prm_the[i]);
        y_diff_the.push_back(dd);
        dd=dd-yn_diff_exp[i]; 
        dd=fabs(dd);
        diff_expthe.push_back(dd);
    } 
    double result;

     
    vector<double> xout_o3; 
    vector<double> yexpout_o3,ydiffout_o3; 
    vector<double> yintout_o3,yprmout_o3; 
    if(myinput.OBJ==3)
    { 
      
        result=0;
        for(int i=0;i<x_diff_exp.size();i++)
        {
            if(x_diff_exp[i]>myinput.OBJ_aout&&x_diff_exp[i]<myinput.OBJ_bout)
                {
                xout_o3.push_back(x_diff_exp[i]);
                yexpout_o3.push_back(y_diff_exp[i]);
                } 
        } 
        interp1d(x_prm_the,y_prm_the,xout_o3,yprmout_o3);
        interp1d(x_xmu,y_xmu,xout_o3,yintout_o3);
         
        int nfit=0;
        for(unsigned i=0;i<=xout_o3.size()-1;++i)
        {
            double yi=0;
            yi=ex*(yintout_o3[i]-yprmout_o3[i]);
            ydiffout_o3.push_back(yi);
        } 
        for(unsigned i=0;i<=xout_o3.size()-1;++i)
        {
            if(xout_o3[i]>myinput.OBJ_afit&&xout_o3[i]<myinput.OBJ_bfit)
            {
            double dd = 0;
             
            dd=ydiffout_o3[i]-yexpout_o3[i];
            dd=dd*dd;
            result=result+dd;
            nfit=nfit+1;
            }
        } 
        result=result*(double)myinput.npar/(double)nfit;
        double OBJ_eps=1;
        double eps_2=OBJ_eps*OBJ_eps;
        result=result/eps_2;
    } 


    if(myinput.OBJ==2){
        result=0;
        for(unsigned i=0;i<=xn.size()-1;++i){
            double dd = 0;
            dd=ex*(yn_xmu[i]-yn_prm_the[i])-yn_diff_exp[i]; 
            dd=dd*dd;
            result=result+dd;
        } 
        result=result*(double)myinput.npar/(double)myinput.ip_num;
        double eps_2=myinput.OBJ_eps*myinput.OBJ_eps;
        result=result/eps_2;
    } 
    if(myinput.OBJ==1){
    result=trapz(diff_expthe,xn[1]-xn[0]);
    double xlen=*(xn.end()-1)-*xn.begin();
    result=result/xlen;} 
     
    if(myinput.OBJ==11)
        {
        result=0;
        double up=0;
        double down=0;
        for(unsigned i=0;i<=xn.size()-1;++i)
            {
            double dd = 0;
            dd=ex*(yn_xmu[i]-yn_prm_the[i])-yn_diff_exp[i]; 
            dd=dd*dd;
            up=up+dd;
            down=down+(yn_diff_exp[i]*yn_diff_exp[i]);
            } 
         result=up/down;
        } 

     
    const char* x_char;
    string title;
    ostringstream ostr;
    ostr<<"saveall/";
     
    ostr<<result<<"struct_";
    for(int i=0;i<=n_fro-1;++i)
        ostr<<x_fro[i]<<"_";
     
    ostr<<"_esdiff_ex_";
    for(int i=0;i<=n-1;++i)
        ostr<<x[i]<<"_";

    string strfix("diffsave.txt");
    title=title+ostr.str()+strfix;
    x_char=title.c_str();
     
     
    ofstream optsave(x_char);
     
    if(myinput.OBJ==1||myinput.OBJ==2||myinput.OBJ==11)
    {
        for(unsigned i=0;i<=xn.size()-1;++i)
            
           optsave<<xn[i]<<"  "<<yn_diff_exp[i]<<"  "<<y_diff_the[i]<<" "<<yn_xmu[i]<<" "<<yn_prm_the[i]<<endl;
    }
    if(myinput.OBJ==3)
    {
        for(unsigned i=0;i<=xout_o3.size()-1;++i)
            optsave<<xout_o3[i]<<"  "<<yexpout_o3[i]<<"  "<<ydiffout_o3[i]<<" "<<yintout_o3[i]<<" "<<yprmout_o3[i]<<endl;
    }

    

    return result;
}



double trapz(vector<double> y,double h){
     
     

    int n=y.size();
    double sum = 0;
    for(int i=1;i<=n-2;++i){
        sum=sum+y[i];
    } 
    sum=y[0]+2*sum+y[n-1];
    sum=sum*(h/2);
    return sum;
} 


} 
#endif  
