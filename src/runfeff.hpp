#ifndef RUNFIT
#define RUNFIT








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
 

#include "subopt.hpp" 
using namespace feff;
#include "subdiff.hpp" 
using namespace feffdiff;
#include "subff2x.hpp"
using namespace ff2x; 
#include "subff2xdiff.hpp"
using namespace ff2xdiff; 

#include "mytool.hpp"
#include "glo.hpp"



 
jmp_buf jb;


using namespace std;




double mycfit(int npar,double* x); 
 
void handler_tiem_out(int xx);

int feffinp(vector<PDB> pdbline);
int feffinp(vector<PDB> pdbline,char workdir_char[CHAR_LEN]);


double mycfit(int npar,double* x){ 
 

 

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
 
    #ifdef MINGW 
    mkdir("feffrun");
    #else
    mkdir("feffrun",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #endif  

strcpy(workdir_char,"feffrun"); 
 
feffinp(pdbline); 
feffinp(pdbline,workdir_char); 
 

 
 
char dirnow_char[CHAR_LEN]; 
getcwd(dirnow_char,sizeof(dirnow_char));
chdir(workdir_char); 

 
#ifndef MINGW
signal(SIGALRM, handler_tiem_out);
alarm(TIME_OUT);  
#endif  
if (setjmp(jb) != 0){

        chdir(dirnow_char); 

         
         
 
         
         
         
         
        result=1e20;
        wrongfile("saveall",x,npar,"timeout");
        return result;

} 
 
 
 
if(myinput.runCore!=0){
system(myinput.feffrun);
}

 

 
chdir(dirnow_char); 
 

char xmupos_char[CHAR_LEN];
memcpy(xmupos_char,workdir_char,CHAR_LEN);
strcat(xmupos_char,"/xmu.dat");
if(access(xmupos_char,0)==-1) {

        char rm_char[CHAR_LEN]=RMDIR; 

        strcat(rm_char,workdir_char);
         
        if(myinput.savedir!=1)system(rm_char);
        result=1e21;
        wrongfile("saveall",x,npar,"noxmu");
        return result;
} 
 
 
if(myinput.runtype==1)result=feff::subopt(npar,x,xmupos_char);

if(myinput.runtype==2)result=feffdiff::subopt(npar,x,xmupos_char);

if(myinput.runtype==5) result=ff2x::subopt(npar,x,xmupos_char);
if(myinput.runtype==6) result=ff2xdiff::subopt(npar,x,xmupos_char);
 
char rm_char[CHAR_LEN]=RMDIR; 
strcat(rm_char,workdir_char);
 
if(myinput.savedir!=1)system(rm_char);
return result;

} 

void handler_tiem_out(int xx){
     
     
     
    longjmp(jb, 1);

} 



int feffinp(vector<PDB> pdbline){
    fstream fm("feffmodel.inp",ios::in);
    string ol;
    vector<string> feffline;
    while(getline(fm,ol)){
        if(ol.length()==0) continue;
        feffline.push_back(ol);
    } 
    ofstream fnew("feff.inp");
    for(vector<string>::iterator i=feffline.begin();i!=feffline.end();++i)
        fnew<<*i<<endl;
    ofstream fnewapp("feff.inp",ios::app);
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i)
        fnewapp<<boost::format("%8.3f   %8.3f   %8.3f %i  %s  %8.3f")%(i->xn)%(i->yn)%(i->zn)%(i->pot)%(i->group)%(i->rn)<<endl;
    fnewapp<<"END"<<endl;

return 0;} 

int feffinp(vector<PDB> pdbline,char workdir_char[CHAR_LEN]){
    char fnew_char[CHAR_LEN];
    fstream fm("feffmodel.inp",ios::in);
    string ol;
    vector<string> feffline;
     
    memcpy(fnew_char,workdir_char,CHAR_LEN);
    strcat(fnew_char,"/feff.inp");
     
    while(getline(fm,ol)){
        if(ol.length()==0) continue;
        feffline.push_back(ol);
    } 
    ofstream fnew(fnew_char);
    for(vector<string>::iterator i=feffline.begin();i!=feffline.end();++i)
        fnew<<*i<<endl;
    ofstream fnewapp(fnew_char,ios::app);
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i)
        fnewapp<<boost::format("%8.3f   %8.3f   %8.3f %i  %s  %8.3f")%(i->xn)%(i->yn)%(i->zn)%(i->pot)%(i->group)%(i->rn)<<endl;
    fnewapp<<"END"<<endl;
return 0;
} 


#endif  
