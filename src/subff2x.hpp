#ifndef SUBFF2X
#define SUBFF2X
namespace ff2x{
#include <iostream>
#include <fstream>
#include<sstream>
#include<vector>

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include <unistd.h> 
 
 
#include "ipmy.hpp"
#include "read_dat.hpp"
#include "nlopt.h"
#include "glo.hpp"

using namespace std;

struct opt_trans{
    string xmupos_string;
    vector<double> x1;
    vector<double> y1;
    vector<double> x2;
    vector<double> y2;
}; 

double trapz(vector<double> y,double h);
double optfun(unsigned n, const double *x, double *grad, void * f_data);
double optfun_save(int n_fro,double* x_fro,int n,const double *x,void * f_data);
int readprm(vector<double>& ,vector<double>&);
int readxmu(vector<double>& ,vector<double>& );
int ff2xinp(double,double);
 
 
 
 
double subopt(int n_fro,double* x_fro,char xmupos_char[]){ 
     

cout<<"in FEFF BOARD"<<endl;
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
myt.xmupos_string=xmupos_char;
 
int npar=4; 
double ra_ori=yn_prm.back()/yn_xmu.back(); 
 
 
double lb[4] = {myinput.lo_es,myinput.lo_ra*ra_ori,myinput.lo_ef_feff,myinput.lo_conv_feff};
double up[4] = {myinput.up_es,myinput.up_ra*ra_ori,myinput.up_ef_feff,myinput.up_conv_feff};
double x0[4] = {myinput.ini_es,myinput.ini_ra*ra_ori,myinput.ini_ef_feff,myinput.ini_conv_feff};   

 
 
 


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
    printf("non-str parameter nlopt failed!\n");
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
    ff2xinp(x[2],x[3]);
    chdir("./feffrun");
    system("rdinp");
    system("ff2x");
    chdir("../");
    opt_trans* tt=(opt_trans *) f_data;
    vector<double> xo_vec,yo_vec;
    readxmu(xo_vec,yo_vec);

 
 
 
 
    vector<double> x0_vec = tt->x1; 
    vector<double> y0_vec = tt->y1;
    string xmupos_string=tt->xmupos_string; 
    const char* xmupos_char=xmupos_string.c_str();

    vector<double> x1_vec ;
    vector<double> y1_vec ;
    vector<double> x2_vec ;
    vector<double> y2_vec ;
    readxmu(x1_vec,y1_vec);
    readprm(x2_vec,y2_vec);
     
    vector<double> xr_vec,yr1_vec,yr2_vec; 
    vector<double> yr0_vec; 
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
     
    if(myinput.OBJ=1){
            
    result=trapz(yr_diff,xr_vec[1]-xr_vec[0]);
    double xlen=*(xr_vec.end()-1)-*xr_vec.begin();
    result=result/xlen;} 


     
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 



    return result;
} 


double optfun_save(int n_fro,double* x_fro,int n,const double *x,void * f_data){
     


    double es = x[0]; 
    double ra = x[1]; 
    double ef=x[2];double cb=x[3];
    ff2xinp(x[2],x[3]);
    chdir("./feffrun");
    system("rdinp");
    system("ff2x");
    chdir("../");
    opt_trans* tt=(opt_trans *) f_data;
    vector<double> xo_vec,yo_vec;
    readxmu(xo_vec,yo_vec);

 
 
 
 
    vector<double> x0_vec = tt->x1; 
    vector<double> y0_vec = tt->y1;
    string xmupos_string=tt->xmupos_string; 
    const char* xmupos_char=xmupos_string.c_str();

    vector<double> x1_vec ;
    vector<double> y1_vec ;
    vector<double> x2_vec ;
    vector<double> y2_vec ;
    readxmu(x1_vec,y1_vec);
    readprm(x2_vec,y2_vec);
     
    vector<double> xr_vec,yr1_vec,yr2_vec; 
    vector<double> yr0_vec; 
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
     
    if(myinput.OBJ=1){
            
    result=trapz(yr_diff,xr_vec[1]-xr_vec[0]);
    double xlen=*(xr_vec.end()-1)-*xr_vec.begin();
    result=result/xlen;} 


     
    const char* x_char;
    string title;
    ostringstream ostr;
    ostr<<"saveall/";
     
    ostr<<result<<"struct_";
    for(int i=0;i<=n_fro-1;++i)
        ostr<<x_fro[i]<<"_";
     
    ostr<<"_es_ratio_";
     
    ostr<<es<<"_"<<ra<<"_"<<ef<<"_"<<cb;

    string strfix("optsave.txt");
    title=title+ostr.str()+strfix;
    x_char=title.c_str();
     
    ofstream optsave(x_char);
    for(unsigned i=0;i<=xr_vec.size()-1;++i)
       optsave<<boost::format("%15.3f %15.3f %15.3f")%(xr_vec[i])%(yr1_vec[i])%(yr2_vec[i])<<endl;
       



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

int readxmu(vector<double>& x_xmu,vector<double>& y_xmu){
     
    vector<vector<double> > dat_xmu;
    read_dat("feffrun/xmu.dat",6,dat_xmu);
    x_xmu=dat_xmu[0];y_xmu=dat_xmu[3];
} 
int readprm(vector<double>& x_prm,vector<double>& y_prm ){
    vector<vector<double> > dat_prm;
    read_dat(myinput.prmname,2,dat_prm);
    x_prm=dat_prm[0];y_prm=dat_prm[1];
} 

int ff2xinp(double ef,double conv){

    fstream fm("./ff2xmodel.inp",ios::in); 
    string ol;
    vector<string> feffline;
     
    while(getline(fm,ol)){
        if(ol.length()==0) continue;
        feffline.push_back(ol);
    } 
    ofstream fnew("./feffrun/feff.inp");
    fnew<<"CORRECTIONS "<<ef<<" "<<conv<<endl;
    for(vector<string>::iterator i=feffline.begin();i!=feffline.end();++i) fnew<<*i<<endl;
     
return 0;
} 


 



} 
#endif  
