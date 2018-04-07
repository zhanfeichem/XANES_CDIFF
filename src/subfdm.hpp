#ifndef SUBFDM
#define SUBFDM

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
 
#define DISsubfdm
#define DEBUGsubfdm

namespace fdm{ 



using namespace std;

struct opt_trans{
    vector<double> x1;
    vector<double> y1;
    vector<double> x2;
    vector<double> y2;
}; 

double trapz(vector<double> y,double h);
int convinp(const double* par_conv);

double convopt();

double convfun(unsigned n, const double *x, double *grad, void * f_data); 
double convfun_save(unsigned n, const double *x, double *grad, void * f_data); 

double subopt();
double subopt_save(); 

double optfun(unsigned n, const double *x, double *grad, void * f_data); 
double optfun_save(int n, const double *x,               void * f_data);

 
 
 
 

double convopt(){


double* lb=(double*)malloc(myinput.n_conv*sizeof(double));
double* up=(double*)malloc(myinput.n_conv*sizeof(double));
double* x0=(double*)malloc(myinput.n_conv*sizeof(double));
 
for(int i=0;i<myinput.n_conv;i++){
    lb[i]=myinput.fdm_low_conv[i];
    up[i]=myinput.fdm_up_conv[i];
    x0[i]=myinput.fdm_x0_conv[i];
}
nlopt_opt opt;
 
 
if(myinput.nstr_ALG=="COBYLA") opt = nlopt_create(NLOPT_LN_COBYLA,myinput.n_conv);
if(myinput.nstr_ALG=="BOBYQA") opt = nlopt_create(NLOPT_LN_BOBYQA,myinput.n_conv);
if(myinput.nstr_ALG=="PRAXIS") opt = nlopt_create(NLOPT_LN_PRAXIS,myinput.n_conv);
if(myinput.nstr_ALG=="NELDERMEAD") opt = nlopt_create(NLOPT_LN_NELDERMEAD,myinput.n_conv);
if(myinput.nstr_ALG=="SBPLX") opt = nlopt_create(NLOPT_LN_SBPLX,myinput.n_conv);
 
if(myinput.nstr_ALG=="ISRES") opt = nlopt_create(NLOPT_GN_ISRES,myinput.n_conv);
if(myinput.nstr_ALG=="ESCH") opt = nlopt_create(NLOPT_GN_ESCH,myinput.n_conv);
if(myinput.nstr_ALG=="DIRECT_L") opt = nlopt_create(NLOPT_GN_DIRECT_L,myinput.n_conv);
if(myinput.nstr_ALG=="CRS2_LM") opt = nlopt_create(NLOPT_GN_CRS2_LM,myinput.n_conv);
 
nlopt_set_lower_bounds(opt,lb);
nlopt_set_upper_bounds(opt,up);
nlopt_set_min_objective(opt,convfun,NULL); 
 
nlopt_set_maxeval(opt,myinput.nopt_conv); 
double minf; /* the minimum objective value, upon return */
 
if (nlopt_optimize(opt, x0, &minf) < 0) {
    printf("conv nlopt failed!\n");
    minf=1e30;
}
else {
     
}
nlopt_destroy(opt);
 
 
double result=0;
result=convfun_save(myinput.n_conv,x0,NULL,NULL);
 
 
return result; 
} 

double convfun(unsigned n, const double *x, double *grad, void * f_data){
#ifdef DIS
if(n!=5)cout<<"The convolution parameter is not 5"<<endl;
#endif  
 
convinp(x);
 
 
 
 
char workdir_char[CHAR_LEN];
strcpy(workdir_char,"fdmnesrun");
char dirnow_char[CHAR_LEN]; 
getcwd(dirnow_char,sizeof(dirnow_char)); 
chdir(workdir_char); 
system(myinput.fdmneswin);
chdir(dirnow_char); 
 
 
double result=0;
result=subopt();
 
 
return result;
 
};

double convfun_save(unsigned n, const double *x, double *grad, void * f_data){
#ifdef DIS
if(n!=5)cout<<"The convolution parameter is not 5"<<endl;
#endif  
 
convinp(x);
 
 
myinput.fdm_parconv=vector<double>(myinput.n_conv); 
for(int i=0;i<myinput.n_conv;i++)myinput.fdm_parconv[i]=x[i];
 
char workdir_char[CHAR_LEN];
strcpy(workdir_char,"fdmnesrun");
char dirnow_char[CHAR_LEN]; 
getcwd(dirnow_char,sizeof(dirnow_char)); 
chdir(workdir_char); 
system(myinput.fdmneswin);
chdir(dirnow_char); 
 
 

double result=0;
result=subopt_save();
 
cout<<"x value is: ";
 
for(int i=0;i<myinput.fdm_ngeom;i++)cout<<myinput.fdm_pargeom[i]<<" ";
cout<<"objective function value: "<<result<<endl;
 
return result;
 
};

double subopt(){ 
     
    vector<vector<double> > dat_xmu;
    vector<vector<double> > dat_prm;
    vector<double> x_xmu,y_xmu,x_prm,y_prm;
    vector<double> xn,yn_xmu,yn_prm;

     
    read_dat_fdm("fdmnesrun/out_conv.txt",2,dat_xmu);
    read_dat(myinput.prmname,2,dat_prm);
    x_xmu=dat_xmu[0];y_xmu=dat_xmu[1];
    x_prm=dat_prm[0];y_prm=dat_prm[1];
 
 
 

 
opt_trans myt; 
myt.x1=x_xmu;
myt.y1=y_xmu;
myt.x2=x_prm;
myt.y2=y_prm;
 
int npar=2;
ip(x_xmu,y_xmu,x_prm,y_prm,xn,yn_xmu,yn_prm,myinput.ip_num); 
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

 
 
 
nlopt_destroy(opt);
 
return minf; 
} 

double subopt_save(){ 
     
    vector<vector<double> > dat_xmu;
    vector<vector<double> > dat_prm;
    vector<double> x_xmu,y_xmu,x_prm,y_prm;
    vector<double> xn,yn_xmu,yn_prm;

     
    read_dat_fdm("fdmnesrun/out_conv.txt",2,dat_xmu);
    read_dat(myinput.prmname,2,dat_prm);
    x_xmu=dat_xmu[0];y_xmu=dat_xmu[1];
    x_prm=dat_prm[0];y_prm=dat_prm[1];
 
 
 

 
opt_trans myt; 
myt.x1=x_xmu;
myt.y1=y_xmu;
myt.x2=x_prm;
myt.y2=y_prm;
 
int npar=2;
ip(x_xmu,y_xmu,x_prm,y_prm,xn,yn_xmu,yn_prm,myinput.ip_num); 
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
 printf("subopt nlopt failed!\n");
}
else {
 
 
 
}

 
optfun_save(2,x0,&myt);
 
nlopt_destroy(opt);
 
return minf; 
} 




int convinp(const double* par_conv){
     

    char workdir_char[CHAR_LEN];
    strcpy(workdir_char,"fdmnesrun"); 
    char fnew_char[CHAR_LEN];
    fstream fm("fdm_conv.inp",ios::in);
    string ol;
    vector<string> feffline;
    memcpy(fnew_char,workdir_char,CHAR_LEN); 
    strcat(fnew_char,"/fdmnes.inp");
    while(getline(fm,ol)){ 
        if(ol.length()==0) continue;
        feffline.push_back(ol);
    } 
    ofstream fnew(fnew_char);
    for(vector<string>::iterator i=feffline.begin();i!=feffline.end();++i)
        fnew<<*i<<endl;
    ofstream fnewapp(fnew_char,ios::app);
     
    fnewapp<<"Efermi"<<endl;
    fnewapp<<boost::format("%3.3f")%(par_conv[0])<<endl;
    fnewapp<<"Gamma_hole"<<endl;
    fnewapp<<boost::format("%3.3f")%(par_conv[1])<<endl;
    fnewapp<<"Gamma_max"<<endl;
    fnewapp<<boost::format("%3.3f")%(par_conv[2])<<endl;
    fnewapp<<"Ecent"<<endl;
    fnewapp<<boost::format("%3.3f")%(par_conv[3])<<endl;
    fnewapp<<"Elarg"<<endl;
    fnewapp<<boost::format("%3.3f")%(par_conv[4])<<endl;
    fnewapp<<"END"<<endl;
return 0;
} 





double optfun_save(int n,const double *x,void * f_data){
     
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
    ip(x1_vec,y1_vec,x2_vec,y2_vec,xr_vec,yr1_vec,yr2_vec,1000);
    for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=fabs(dd);
        yr_diff.push_back(dd);
    } 
    double result;
    result=trapz(yr_diff,xr_vec[1]-xr_vec[0]);
    double xlen=*(xr_vec.end()-1)-*xr_vec.begin();
    result=result/xlen;

     
    const char* x_char;
    string title;
    ostringstream ostr;
    ostr<<"saveall/";
     
    ostr<<result<<"struct_"; 
    for(int i=0;i<=myinput.fdm_ngeom-1;++i)
        ostr<<myinput.fdm_pargeom[i]<<"_";
     
    ostr<<"conv_"; 
    for(int i=0;i<=myinput.n_conv-1;++i)
        ostr<<myinput.fdm_parconv[i]<<"_";
     
    ostr<<"_es_ratio_";
    for(int i=0;i<=n-1;++i)
        ostr<<x[i]<<"_";

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
    ip(x1_vec,y1_vec,x2_vec,y2_vec,xr_vec,yr1_vec,yr2_vec,1000);

    for(unsigned i=0;i<=xr_vec.size()-1;++i){
        double dd = 0;
        dd=yr1_vec[i]-yr2_vec[i];
        dd=fabs(dd);
        yr_diff.push_back(dd);
    } 
    double result;
    result=trapz(yr_diff,xr_vec[1]-xr_vec[0]);
    double xlen=*(xr_vec.end()-1)-*xr_vec.begin();
     
 
 
 
 
 
 

    return result/xlen;
}





} 
#endif  
