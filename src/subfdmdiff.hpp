#ifndef SUBFDMDIFF
#define SUBFDMDIFF

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
 
#define DISsubfdmdiff
#define DEBUGsubfdmdiff

namespace fdmdiff{ 


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


double convopt();
double convfun(unsigned n, const double *x, double *grad, void * f_data); 
double convfun_save(unsigned n, const double *x, double *grad, void * f_data); 
int convinp(const double* par_conv);
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
    printf("nlopt failed!\n");
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
 
 
return subopt();
 
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
    vector<double> x_xmu,y_xmu;
    read_dat_fdm("fdmnesrun/out_conv.txt",2,dat_xmu);
    x_xmu=dat_xmu[0];y_xmu=dat_xmu[1];
     
    vector<vector<double> > dat_prm_the;
    vector<double> x_prm_the,y_prm_the;
    read_dat(myinput.prm_the_char,2,dat_prm_the);

    x_prm_the=dat_prm_the[0];y_prm_the=dat_prm_the[1];
     
    vector<vector<double> > dat_diff_exp; 
    vector<double> x_diff_exp,y_diff_exp;
    read_dat(myinput.diff_exp_char,2,dat_diff_exp);
    x_diff_exp=dat_diff_exp[0];y_diff_exp=dat_diff_exp[1];
 
 
 
 

 
opt_trans myt; 
myt.x_xmu=x_xmu;
myt.y_xmu=y_xmu;
myt.x_prm_the=x_prm_the;
myt.y_prm_the=y_prm_the;
myt.x_diff_exp=x_diff_exp;
myt.y_diff_exp=y_diff_exp;
myt.es_prm=myinput.es_prm;
myt.ra_prm=myinput.ra_prm;
 
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
 
nlopt_destroy(opt);
return minf; 
}

double subopt_save(){ 
     
     
    vector<vector<double> > dat_xmu;
    vector<double> x_xmu,y_xmu;
    read_dat_fdm("fdmnesrun/out_conv.txt",2,dat_xmu);
    x_xmu=dat_xmu[0];y_xmu=dat_xmu[1];
     
    vector<vector<double> > dat_prm_the;
    vector<double> x_prm_the,y_prm_the;
    read_dat(myinput.prm_the_char,2,dat_prm_the);
    x_prm_the=dat_prm_the[0];y_prm_the=dat_prm_the[1];
     
    vector<vector<double> > dat_diff_exp; 
    vector<double> x_diff_exp,y_diff_exp;
    read_dat(myinput.diff_exp_char,2,dat_diff_exp);
    x_diff_exp=dat_diff_exp[0];y_diff_exp=dat_diff_exp[1];

 
 
 

 
opt_trans myt; 
myt.x_xmu=x_xmu;
myt.y_xmu=y_xmu;
myt.x_prm_the=x_prm_the;
myt.y_prm_the=y_prm_the;
myt.x_diff_exp=x_diff_exp;
myt.y_diff_exp=y_diff_exp;
myt.es_prm=myinput.es_prm;
myt.ra_prm=myinput.ra_prm;
 

 
 
 
 

 
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

optfun_save(2,x0,&myt); 

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


     
    #ifdef DEBUGsubfdmdiff
    ofstream oo1("./debug/diff1.txt");
    for(int i=0;i<xn.size();i++) oo1<<xn[i]<<" "<<yn_xmu[i]<<" "<<yn_diff_exp[i]<<" "<<yn_prm_the[i]<<endl;
    #endif  



    vector<double> diff_expthe; 
    diff_expthe.clear();
    for(unsigned i=0;i<=xn.size()-1;++i){
        double dd = 0;
        dd=ex*(yn_xmu[i]-yn_prm_the[i])-yn_diff_exp[i]; 

        dd=fabs(dd);
        diff_expthe.push_back(dd);
    } 
    double result;
    result=trapz(diff_expthe,xn[1]-xn[0]);
    double xlen=*(xn.end()-1)-*xn.begin();
    result=result/xlen;
    return result;
}
double optfun_save(int n,const double *x,void * f_data){
     
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
    result=trapz(diff_expthe,xn[1]-xn[0]);
    double xlen=*(xn.end()-1)-*xn.begin();
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
     
    ostr<<"_esdiff_ex_";
    for(int i=0;i<=n-1;++i)
        ostr<<x[i]<<"_";

    string strfix("diffsave.txt");
    title=title+ostr.str()+strfix;
    x_char=title.c_str();
     
    ofstream optsave(x_char);
    for(unsigned i=0;i<=xn.size()-1;++i)
       optsave<<xn[i]<<"  "<<y_diff_the[i]<<"  "<<yn_diff_exp[i]<<" "<<yn_xmu[i]<<" "<<yn_prm_the[i]<<endl;
        
    

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



















} 
#endif  
