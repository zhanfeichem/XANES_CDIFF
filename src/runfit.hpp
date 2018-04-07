#ifndef RUNFIT
#define RUNFIT




namespace feffmain{



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


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


int n_change=1;
change_shell(n_change,'B',x[0],pdbline);
change_shell(n_change,'C',x[1],pdbline);
change_shell(n_change,'D',x[1],pdbline);
change_shell(n_change,'E',x[2],pdbline);
change_shell(n_change,'F',x[2],pdbline);






 
save_pdb(pdbline,npar,x);
 
char workdir_char[CHAR_LEN];
myworkdir(workdir_char); 
 
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
 
 

 

 

 
chdir(dirnow_char); 
 

char xmupos_char[CHAR_LEN];
memcpy(xmupos_char,workdir_char,CHAR_LEN);
strcat(xmupos_char,"/xmu.dat");
if(access(xmupos_char,0)==-1) {

        char rm_char[CHAR_LEN]=RMDIR; 

        strcat(rm_char,workdir_char);

        system(rm_char);
        result=1e21;
        wrongfile("saveall",x,npar,"noxmu");
        return result;
} 
 
 
result=feffdiff::subopt(npar,x,xmupos_char);
 
char rm_char[CHAR_LEN]=RMDIR; 
strcat(rm_char,workdir_char);

 
 

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



}

#endif  
