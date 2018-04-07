#include <iostream>
#include <fstream>
#include <string>
 
#include "runfdm.hpp" 
#include "runfeff.hpp" 
#include "dl_fdm.hpp" 
 
 
 
 
 
#include "read_dat.hpp"
#include <cstring>
 
#ifdef USENOMAD
#include "nomad.hpp"
#endif
using namespace std;
 
struct trans_npar{
int npar;
}; 

int run_nlopt();
int run_nomadlib(int argc , char ** argv);
double myfunc(unsigned , const double *, double *, void *);
double fdm_myfunc(unsigned n, const double *x, double *grad, void *my_func_data);
int run_dl(); 

int main ( int argc , char ** argv ) {
 
inputfile();
cout<<"after input"<<endl;
if(myinput.OPT=="NOMAD"){
int npar = myinput.npar; 
cout<<"npar is"<<npar;
double f = 1e20;
double x[npar];
if ( argc >= 2 ) {
ifstream in ( argv[1] ); 
for ( int i = 0 ; i <=npar-1 ; i++ ) in >> x[i]; 
 
if(myinput.runtype==1 or myinput.runtype==2)f = mycfit(npar,x); 

if (in.fail()) f = 1e19;
in.close();
}
cout << f << endl; 
return 0;} 

if(myinput.OPT=="NLOPT") run_nlopt(); 
#ifdef USENOMAD
if(myinput.OPT=="NOMADLIB") run_nomadlib(argc,argv);
#endif  

if(myinput.OPT=="DL") run_dl();

} 

 
int run_dl(){
 

double d1[11]={-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5};
double d2[11]={-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5};
double d3[11]={-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5};
int ngrid=11;
double f=0.0;
int npar=3;
double x[3]={0,0,0};
for(int i=0;i<ngrid;++i)
    for(int j=0;j<ngrid;++j)
        for(int k=0;k<ngrid;++k){
            x[0]=d1[i];x[1]=d2[j];x[2]=d3[k];

            dl_mycfit_fdm(npar,x); 
        } 

 
 
 
 
 

}


 
int run_nlopt (  ) {


 
int npar = myinput.npar; 
int nval=myinput.NLOPT_nval; 
string OPT=myinput.NLOPT_ALG;
double lb[npar];
for(int i=0;i<npar;i++)lb[i]=myinput.NLOPT_xl[i];
double ub[npar];
for(int i=0;i<npar;i++)ub[i]=myinput.NLOPT_xu[i];
double x[npar]; 
for(int i=0;i<npar;i++)x[i]=myinput.NLOPT_x0[i];
 
nlopt_opt opt;
 
if(OPT=="DIRECT_L") opt = nlopt_create(NLOPT_GN_DIRECT_L,npar);
if(OPT=="CRS2_LM") opt = nlopt_create(NLOPT_GN_CRS2_LM,npar);
if(OPT=="ISRES") opt = nlopt_create(NLOPT_GN_ISRES,npar);
if(OPT=="ESCH") opt = nlopt_create(NLOPT_GN_ESCH,npar);
 
if(OPT=="NELDERMEAD") opt = nlopt_create(NLOPT_LN_NELDERMEAD,npar);
if(OPT=="SBPLX") opt = nlopt_create(NLOPT_LN_SBPLX,npar);
if(OPT=="PRAXIS") opt = nlopt_create(NLOPT_LN_PRAXIS,npar);
if(OPT=="BOBYQA") opt = nlopt_create(NLOPT_LN_BOBYQA,npar);
if(OPT=="NEWUOA_BOUND") opt = nlopt_create(NLOPT_LN_NEWUOA_BOUND,npar);
if(OPT=="COBYLA") opt = nlopt_create(NLOPT_LN_COBYLA,npar);

 
nlopt_set_lower_bounds(opt, lb);
nlopt_set_upper_bounds(opt,ub);
nlopt_set_maxeval(opt,nval); 
nlopt_set_xtol_rel(opt,1e-10); 
 
 
 
trans_npar tnpar;
tnpar.npar=myinput.npar;
if(myinput.runtype==1||myinput.runtype==2||myinput.runtype==5||myinput.runtype==6)nlopt_set_min_objective(opt,myfunc,&tnpar);
if(myinput.runtype==3||myinput.runtype==4)nlopt_set_min_objective(opt,fdm_myfunc,&tnpar);
 
double minf;  

if (nlopt_optimize(opt, x, &minf) < 0) {printf(" structure search nlopt failed!\n");}
else { printf("found minimum at f = %0.10g\n", minf);}
 
cout<<"The x value is ";
for(int i;i<myinput.npar;i++)cout<<x[i]<<" ";
cout<<endl;
}

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data){
 
trans_npar* tnpar=(trans_npar*)my_func_data; 
int npar=tnpar->npar;
double xin[npar]; 
for(int i=0;i<npar;i++)xin[i]=x[i];
return mycfit(npar,xin);
}

double fdm_myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
 
trans_npar* tnpar=(trans_npar*)my_func_data; 
int npar=tnpar->npar;
double xin[npar]; 
for(int i=0;i<npar;i++)xin[i]=x[i];
return mycfit_fdm(npar,xin);

}

#ifdef USENOMAD
 
class My_Evaluator : public NOMAD::Evaluator {
public:
  My_Evaluator  ( const NOMAD::Parameters & p ) :
    NOMAD::Evaluator ( p ) {}

  ~My_Evaluator ( void ) {}

bool eval_x ( NOMAD::Eval_Point   & x, const NOMAD::Double & h_max,bool& count_eval   ) const{

 
 
int npar=myinput.npar;
 
double xd[npar];
for(int i=0;i<npar;i++) xd[i]=x[i].value();
 
double result=mycfit(npar,xd);
 
NOMAD::Double result_nomad=result; 
 
x.set_bb_output  ( 0 , result_nomad );
 
count_eval = true;  
return true;        
} 

}; 

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int run_nomadlib ( int argc , char ** argv) {

   
  NOMAD::Display out ( std::cout );
  out.precision ( NOMAD::DISPLAY_PRECISION_STD );

  try {

     
    NOMAD::begin ( argc , argv );
     
    NOMAD::Parameters p ( out );
     
    int npar=myinput.npar;
    p.set_DIMENSION (npar);              
     
    vector<NOMAD::bb_output_type> bbot (1);  
    bbot[0] = NOMAD::OBJ;                    
 
 
    p.set_BB_OUTPUT_TYPE ( bbot );

 
    p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );
 
NOMAD::Point x0(npar);
NOMAD::Point low(npar);
NOMAD::Point ub(npar);
 
for(int i=0;i<npar;++i) x0[i]=myinput.NLOPT_x0[i];
for(int i=0;i<npar;++i) low[i]=myinput.NLOPT_xl[i];
for(int i=0;i<npar;++i) ub[i]=myinput.NLOPT_xu[i];





p.set_X0(x0);   
p.set_LOWER_BOUND (low);  
p.set_UPPER_BOUND (ub);

    p.set_MAX_BB_EVAL (myinput.NLOPT_nval);      
    p.set_DISPLAY_DEGREE(2);
    p.set_SOLUTION_FILE("sol.txt");

     
    p.check();

     
    My_Evaluator ev   ( p );

     
    NOMAD::Mads mads ( p , &ev );
    mads.run();
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  NOMAD::Slave::stop_slaves ( out );
  NOMAD::end();
  return EXIT_SUCCESS;
} 

#endif  


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 




 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 


