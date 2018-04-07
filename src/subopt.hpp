#ifndef SUBOPT
#define SUBOPT
namespace feff{
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

using namespace std;

struct opt_trans{
    vector<double> x1;
    vector<double> y1;
    vector<double> x2;
    vector<double> y2;
}; 

double trapz(vector<double> y,double h);
double optfun(unsigned n, const double *x, double *grad, void * f_data);
double optfun_save(int n_fro,double* x_fro,int n,const double *x,void * f_data);
 
 
 
 
double subopt(int n_fro,double* x_fro,char xmupos_char[]){ 
     


    ofstream before("before.txt");
    vector<vector<double> > dat_xmu;
    vector<vector<double> > dat_prm;
    vector<double> x_xmu,y_xmu,x_prm,y_prm;
    vector<double> xn,yn_xmu,yn_prm;

     
    read_dat(xmupos_char,6,dat_xmu);
    read_dat(myinput.prmname,2,dat_prm);

    x_xmu=dat_xmu[0];y_xmu=dat_xmu[3];
    x_prm=dat_prm[0];y_prm=dat_prm[1];
     
 
 
 

    ip(x_xmu,y_xmu,x_prm,y_prm,xn,yn_xmu,yn_prm,myinput.ip_num);
    for(unsigned i=0;i<=xn.size()-1;++i)
        before<<xn[i]<<"  "<<yn_xmu[i]<<"  "<<yn_prm[i]<<endl;
 
opt_trans myt; 
myt.x1=x_xmu;
myt.y1=y_xmu;
myt.x2=x_prm;
myt.y2=y_prm;
 
int npar=2;
double ra_ori=yn_prm.back()/yn_xmu.back(); 
 
 
double lb[2] = {myinput.lo_es,myinput.lo_ra*ra_ori};
double up[2] = {myinput.up_es,myinput.up_ra*ra_ori};
double x0[2] = {myinput.ini_es,myinput.ini_ra*ra_ori};   
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
     
    double es = x[0]; 
    double ra = x[1]; 
    opt_trans* tt=(opt_trans *) f_data;
    vector<double> aa(5,5);
    vector<double> x1_vec = tt->x1; 
    vector<double> y1_vec = tt->y1;
    vector<double> x2_vec = tt->x2; 
    vector<double> y2_vec = tt->y2;
    vector<double> xr_vec,yr1_vec,yr2_vec; 
    vector<double> yr_diff;
     
    for(vector<double>::iterator i=x1_vec.begin();i!=x1_vec.end();++i)
        *i=*i+es;
    for(vector<double>::iterator i=y1_vec.begin();i!=y1_vec.end();++i)
        *i=(*i)*ra;
     
    ip(x1_vec,y1_vec,x2_vec,y2_vec,xr_vec,yr1_vec,yr2_vec,myinput.ip_num);

    for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=fabs(dd);
        yr_diff.push_back(dd);
    } 
    double result; 
     
    vector<double> xexp_o3=x2_vec; 
         
    vector<double> yexp_o3=y2_vec; 

    vector<double> xfeff_o3=x1_vec; 
    vector<double> yfeff_o3=y1_vec;
     
    vector<double> xout_o3;
    vector<double> yexpout_o3;
    vector<double> ytheout_o3; 
    if(myinput.OBJ==3){

        for(int i=0;i<xexp_o3.size();i++){
                if(xexp_o3[i]>myinput.OBJ_aout&&xexp_o3[i]<myinput.OBJ_bout){
                    xout_o3.push_back(xexp_o3[i]);
                    yexpout_o3.push_back(yexp_o3[i]);
                }

        } 
         
        interp1d(xfeff_o3,yfeff_o3,xout_o3,ytheout_o3); 
         
        int nfit=0; 
        for(unsigned i=0;i<=xout_o3.size()-1;++i){
            if(xout_o3[i]>myinput.OBJ_afit&&xout_o3[i]<myinput.OBJ_bfit){ 
                double dd = 0;
                dd=ytheout_o3[i]-yexpout_o3[i];
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
        for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=dd*dd;
        result=result+dd;
        } 

    result=result*(double)myinput.npar/(double)myinput.ip_num;
    double eps_2=myinput.OBJ_eps*myinput.OBJ_eps;
     result=result/eps_2;
    } 
     
    if(myinput.OBJ==1){
    result=trapz(yr_diff,xr_vec[1]-xr_vec[0]);
    double xlen=*(xr_vec.end()-1)-*xr_vec.begin();
    result=result/xlen;
    }
     
    if(myinput.OBJ==11){

         
        result=0;
        double up=0;
        double down=0;
        for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=dd*dd;
        double dd2=0;
        dd2=yr2_vec[i]*yr2_vec[i];
        up=up+dd;
        down=down+dd2;
        } 
        result=up/down;
     } 


    return result;
}



double optfun_save(int n_fro,double* x_fro,int n,const double *x,void * f_data){
     
     
    double es = x[0]; 
    double ra = x[1]; 
    opt_trans* tt=(opt_trans *) f_data;
    vector<double> x1_vec = tt->x1;
    vector<double> y1_vec = tt->y1;
    vector<double> x2_vec = tt->x2;
    vector<double> y2_vec = tt->y2;
    vector<double> xr_vec,yr1_vec,yr2_vec; 
    vector<double> yr_diff;
     
    for(vector<double>::iterator i=x1_vec.begin();i!=x1_vec.end();++i)
        *i=*i+es;
    for(vector<double>::iterator i=y1_vec.begin();i!=y1_vec.end();++i)
        *i=(*i)*ra;
    ip(x1_vec,y1_vec,x2_vec,y2_vec,xr_vec,yr1_vec,yr2_vec,myinput.ip_num);
    for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=fabs(dd);
        yr_diff.push_back(dd);
    } 

    double result;
     
    vector<double> xexp_o3=x2_vec; 
         
    vector<double> yexp_o3=y2_vec; 

    vector<double> xfeff_o3=x1_vec; 
    vector<double> yfeff_o3=y1_vec;
     
    vector<double> xout_o3;
    vector<double> yexpout_o3;
    vector<double> ytheout_o3; 
       if(myinput.OBJ==3){

        for(int i=0;i<xexp_o3.size();i++){
                if(xexp_o3[i]>myinput.OBJ_aout&&xexp_o3[i]<myinput.OBJ_bout){
                    xout_o3.push_back(xexp_o3[i]);
                    yexpout_o3.push_back(yexp_o3[i]);
                }

        } 
         
        interp1d(xfeff_o3,yfeff_o3,xout_o3,ytheout_o3); 
         
        int nfit=0; 
        for(unsigned i=0;i<=xout_o3.size()-1;++i){
            if(xout_o3[i]>myinput.OBJ_afit&&xout_o3[i]<myinput.OBJ_bfit){ 
                double dd = 0;
                dd=ytheout_o3[i]-yexpout_o3[i];
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
        for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=dd*dd;
        result=result+dd;
        } 
     result=result*(double)myinput.npar/(double)myinput.ip_num;
     double eps_2=myinput.OBJ_eps*myinput.OBJ_eps;
     result=result/eps_2;
     } 
     
    if(myinput.OBJ==1){
            
    result=trapz(yr_diff,xr_vec[1]-xr_vec[0]);
    double xlen=*(xr_vec.end()-1)-*xr_vec.begin();
    result=result/xlen;} 
     
    if(myinput.OBJ==11){

         
        result=0;
        double up=0;
        double down=0;
        for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=dd*dd;
        double dd2=0;
        dd2=yr2_vec[i]*yr2_vec[i];
        up=up+dd;
        down=down+dd2;
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
     
    ostr<<"_es_ratio_";
    for(int i=0;i<=n-1;++i)
        ostr<<x[i]<<"_";

    string strfix("optsave.txt");
    title=title+ostr.str()+strfix;
    x_char=title.c_str();
     
    ofstream optsave(x_char);
    if(myinput.OBJ==1||myinput.OBJ==2||myinput.OBJ==11)
    {
         
         

    for(unsigned i=0;i<=xr_vec.size()-1;++i)
        
       optsave<<boost::format("%15.3f %15.3f %15.3f")%(xr_vec[i])%(yr1_vec[i])%(yr2_vec[i])<<endl;
    } 

    if(myinput.OBJ==3){
             
        for(unsigned i=0;i<=xout_o3.size()-1;++i)
        optsave<<boost::format("%15.3f %15.3f %15.3f")%(xout_o3[i])%(yexpout_o3[i])%(ytheout_o3[i])<<endl;

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
