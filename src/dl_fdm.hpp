#ifndef DL_FDM
#define DL_FDM
using namespace std;
  
#include <iostream>
#include <fstream>
#include<sstream>
#include<vector>
#include<map>
#include <string>
  
#include<boost/tokenizer.hpp>  
  
#include<boost/format.hpp>  
  
#include<stdlib.h>  
#include<math.h>  
#include <string.h>  
#include<time.h>   
#include<stdlib.h>
#include <unistd.h>  
#include <setjmp.h>
#include<signal.h>
#include <sys/stat.h>  
#include <sys/types.h>  
  
  


#include "subfdmdiff.hpp"  
using namespace fdmdiff;
#include "subfdm.hpp"
using namespace fdm;
  
#include "mytool.hpp"
#include "glo.hpp"
#include "read_dat.hpp"  
double dl_mycfit_fdm(int ,double* );
  
  
double dl_save(int npar,double* x);

double dl_mycfit_fdm(int npar,double* x){
  
double result = 1e21;
vector<PDB> pdbline;
getpdb(pdbline);

  
change_shell_ini(pdbline);

  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
if(myinput.geom_model==1){
int n_change=myinput.geom_ngroup;

for(int i=0;i<myinput.geom_n;i++){
        int nn=myinput.geom_xi[i];
        double change[3]={x[nn],0,0};
        double atom[3]={myinput.geom_atom[i][0],myinput.geom_atom[i][1],myinput.geom_atom[i][2]};
        change_lig(n_change,myinput.geom_group[i],change,atom,pdbline);
}  

}  
else if(myinput.geom_model==11){
int n_change=myinput.geom_ngroup;

for(int i=0;i<myinput.geom_n;i++){

        int nn=myinput.geom_xi[i];

        double change[3]={x[nn],0,0};

          
        double atom[3]={myinput.geom_atom[i][0],myinput.geom_atom[i][1],myinput.geom_atom[i][2]};
        change_lig_val(n_change,myinput.geom_group[i],change,atom,pdbline);
}  

}  

else if(myinput.geom_model==2){
int n_change=myinput.geom_ngroup;
for(int i=0;i<myinput.geom_n;i++){
        int ri=myinput.geom_r[i];
        int thetai=myinput.geom_theta[i];
        int phii=myinput.geom_phi[i];
          
        double change[3]={x[ri],x[thetai],x[phii]};
          
          
        double atom[3]={myinput.geom_atom[i][0],myinput.geom_atom[i][1],myinput.geom_atom[i][2]};
        change_lig_val(n_change,myinput.geom_group[i],change,atom,pdbline);
}  

}  





  
save_pdb(pdbline,npar,x);
  
char workdir_char[CHAR_LEN];
  
  
strcpy(workdir_char,"fdmnesrun");
  
  
fdmnesinp(pdbline,workdir_char);  
  

  
  
char dirnow_char[CHAR_LEN];  
getcwd(dirnow_char,sizeof(dirnow_char));
chdir(workdir_char);  

  

if(myinput.runCore==1)system(myinput.fdmneswin);


  
chdir(dirnow_char);  
  
char xmupos_char[CHAR_LEN];
memcpy(xmupos_char,workdir_char,CHAR_LEN);
strcat(xmupos_char,"/out.txt");  
if(access(xmupos_char,0)==-1) {

  
  
  
  
  
        result=1e22;
        wrongfile("saveall",x,npar,"noxmu");
        remove("./fdmnesrun/out.bk");
        remove("./fdmnesrun/out_conv.bk");
        return result;
}  

  
myinput.fdm_ngeom=npar;
  
myinput.fdm_pargeom=vector<double>(npar);
for(int i=0;i<npar;i++)myinput.fdm_pargeom[i]=x[i];

  
dl_save(npar,x);
  
remove("./fdmnesrun/out.bk");
remove("./fdmnesrun/out_conv.bk");
rename("./fdmnesrun/out.txt","./fdmnesrun/out.bk");
rename("./fdmnesrun/out_conv.txt","./fdmnesrun/out_conv.bk");
  
return 0;


}  




double dl_save(int npar,double* x){
    vector<vector<double> > dat_xmu;
    vector<double> x_xmu,y_xmu;
      
    read_dat_fdm("fdmnesrun/out_conv.txt",2,dat_xmu);  
    x_xmu=dat_xmu[0];y_xmu=dat_xmu[1];

      
    ofstream dlsave("dl_spc.txt",ios::app);
    for(unsigned i=0;i<=y_xmu.size()-1;++i)
       dlsave<<boost::format("%15.8f")%(y_xmu[i]);
    dlsave<<endl;  
      
    ofstream dlpar("dl_par.txt",ios::app);
    for(unsigned i=0;i<npar;i++){
        dlpar<<x[i]<<" ";
        std::cout<<x[i]<<"  ";
    }
    dlpar<<endl;
    cout<<endl;

}













#endif   
