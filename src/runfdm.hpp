#ifndef RUNFDM
#define RUNFDM
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
double mycfit_fdm(int npar,double* x);
 
int fdmnesinp(vector<PDB> pdbline,char workdir_char[CHAR_LEN]);

double mycfit_fdm(int npar,double* x){ 
 
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


else if(myinput.geom_model==3){

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
        return result;
} 

 
myinput.fdm_ngeom=npar;
 
myinput.fdm_pargeom=vector<double>(npar);
for(int i=0;i<npar;i++)myinput.fdm_pargeom[i]=x[i];

 
 
 
 
if(myinput.runtype==3)result=fdm::convopt();
if(myinput.runtype==4)result=fdmdiff::convopt();
 
remove("./fdmnesrun/out.bk");
rename("./fdmnesrun/out.txt","./fdmnesrun/out.bk");
return result;



} 


int fdmnesinp(vector<PDB> pdbline,char workdir_char[CHAR_LEN]){
    char fnew_char[CHAR_LEN];
    fstream fm("fdmnesmodel.inp",ios::in);
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
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i)
        fnewapp<<boost::format("%i %8.3f   %8.3f   %8.3f !%s %8.3f")%(i->pot)%(i->xn)%(i->yn)%(i->zn)%(i->group)%(i->rn)<<endl;
         
    fnewapp<<"END"<<endl;
return 0;
} 















#endif  
