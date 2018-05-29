#ifndef GLO
#define GLO

  
  
#define MINGW
  
  

#ifdef CYGWIN
#include<io.h>
#endif   
#ifdef CENTOS
#include <sys/io.h>  
#endif   

  
#ifdef MINGW  
#define RMDIR "rd /s/q "  
  
#include<io.h>
#else
#define RMDIR "rm -rf "  
  
  
#endif   
  
#include<string.h>
  
#include<map>
#include<vector>
#include<fstream>
using namespace std;
  
#include<boost/tokenizer.hpp>  
  
#include<boost/format.hpp>  
  
#include<cstdlib>

#define TIME_OUT 1200  
#define CHAR_LEN 50  
  

  
  
  
  
  
  


#define NLOPT_CONV NLOPT_LN_BOBYQA
#define NLOPT_SUB  NLOPT_LN_COBYLA
  
  
  
  
#define GEOM_N 20
  
struct par_str{
    int ngeom;
    double* pargeom;
    int nconv;
    double* parconv;
    double* parsub;
};  
struct Myin{
      
    int display;
    int runtype;
    int npar;
    int runCore;  
    int savedir;
      
    int OBJ;  
    double OBJ_eps;  
    double OBJ_afit,OBJ_bfit;  
    double OBJ_aout,OBJ_bout;  
      
    string OPT;
    string NLOPT_ALG;
    int NLOPT_nval;
    double NLOPT_etol;
    vector<double> NLOPT_x0;
    vector<double> NLOPT_xl;
    vector<double> NLOPT_xu;
      
    char indir[CHAR_LEN];  

    char feffbat[CHAR_LEN];
    char feffrun[CHAR_LEN];
    char fdmneswin[CHAR_LEN];
    char pdbname[CHAR_LEN];
    char prmname[CHAR_LEN];  
      
    char prm_the_char[CHAR_LEN];  
    char diff_exp_char[CHAR_LEN];  

    int  ip_num;  
    string nstr_ALG;  
    int  nstr_num;  

    double up_es,lo_es;  
    double up_ra,lo_ra;
    double ini_es,ini_ra;
    double ini_ef_feff,lo_ef_feff,up_ef_feff;  
    double ini_conv_feff,lo_conv_feff,up_conv_feff;  
      
    double es_prm;  
    double ra_prm;  
    double ef_feff_prm,conv_feff_prm;  
    double up_es_diff,lo_es_diff;  
    double up_ex_diff,lo_ex_diff;  
    double ini_es_diff,ini_ex_diff;
    double** atom_array;
    map<string,int> potmap;
      
    int nopt_conv;  
    int n_conv;  

  
  
  
    vector<double> fdm_up_conv;
    vector<double> fdm_low_conv;
    vector<double> fdm_x0_conv;
      
    vector<double> fdm_parconv;
    int fdm_ngeom;
    vector<double> fdm_pargeom;


  
par_str mypar;

  
int geom_model;  
int geom_ngroup;  
int geom_n;  
double geom_atom[GEOM_N][3];  
char  geom_group[GEOM_N];
  
int  geom_xi[GEOM_N];
  
int geom_r[GEOM_N];
int geom_theta[GEOM_N];
int geom_phi[GEOM_N];
};  
  
  
Myin myinput;
par_str mypar;
int constr(char* s1,char* s2,char* out){
      
    memcpy(out,s1,CHAR_LEN);
    strcat(out,s2);
return 0;
}  

int readinp(const char* ,vector<vector<string> >& );
int dis_vec2d(vector<vector<string> > );  
  
int getkey_str(const string , vector<vector<string> >& ,string& );
int getkey_double(const string , vector<vector<string> >& ,double& );
int getkey_int(const string,  vector<vector<string> >& ,int& );
  
int getkey_double3(const string ,  vector<vector<string> >& ,double& ,double & ,double & );
int get_vector_double(int ,const string , vector<vector<string> >& ,vector<double>& );
  
int getpot(int, vector<vector<string> >& ,map<string,int>& );
  
int getatom(int,vector<vector<string> >& ,vector<vector<double> >& );
int get_geom_group(int , vector<vector<string> >& ,vector<char>& );
int get_geom_xi(int , vector<vector<string> >&,vector<int>& );
int get_geom_xngro(const string ,int , vector<vector<string> >& ,vector<int>& );
  
int input_general(vector<vector<string> >&);  
int input_nlopt(vector<vector<string> >&);
int input_feffdiff(vector<vector<string> >&);
int input_feffprm(vector<vector<string> >&);
int input_ff2x(vector<vector<string> >&);
int input_fdmprm(vector<vector<string> >&);


int input_geom_bond(vector<vector<string> >&);
int input_geom_3d(vector<vector<string> >&);
  
int inputfile(){
    string strin;

    const char* finp="input.inp";
    vector<vector<string> > input_vec2d;
    readinp(finp,input_vec2d);
    dis_vec2d(input_vec2d);  

      
    myinput.display=0;  
    getkey_int("display",input_vec2d,myinput.display);
      
    getkey_int("npar",input_vec2d,myinput.npar);
    myinput.runCore=1;  
    getkey_int("runCore",input_vec2d,myinput.runCore);
    myinput.savedir=0;  
    getkey_int("savedir",input_vec2d,myinput.savedir);


      
    myinput.OPT="NOMAD";  
    getkey_str("OPT",input_vec2d,myinput.OPT);
    if(myinput.OPT=="NLOPT"||myinput.OPT=="NOMADLIB") input_nlopt(input_vec2d);  

      
      
    getkey_int("OBJ",input_vec2d,myinput.OBJ);
    if(myinput.OBJ==2)getkey_double("OBJ_eps",input_vec2d,myinput.OBJ_eps);  
    if(myinput.OBJ==3)getkey_double("OBJ_afit",input_vec2d,myinput.OBJ_afit);  
    if(myinput.OBJ==3)getkey_double("OBJ_bfit",input_vec2d,myinput.OBJ_bfit);  
    if(myinput.OBJ==3)getkey_double("OBJ_aout",input_vec2d,myinput.OBJ_aout);  
    if(myinput.OBJ==3)getkey_double("OBJ_bout",input_vec2d,myinput.OBJ_bout);

      

    getkey_int("runtype",input_vec2d,myinput.runtype);
    if(myinput.runtype==1||myinput.runtype==5) input_feffprm(input_vec2d);  
    if(myinput.runtype==2||myinput.runtype==6) input_feffdiff(input_vec2d);  
    if(myinput.runtype==5||myinput.runtype==6) input_ff2x(input_vec2d);  
      

    if(myinput.runtype==3) input_feffprm(input_vec2d);  
    if(myinput.runtype==3) input_fdmprm(input_vec2d);  
    if(myinput.runtype==4) input_feffdiff(input_vec2d);
    if(myinput.runtype==4) input_fdmprm(input_vec2d);  

    getkey_int("geom_model",input_vec2d,myinput.geom_model);
    if(myinput.geom_model==1||myinput.geom_model==11) input_geom_bond(input_vec2d);
    if(myinput.geom_model==2) input_geom_3d(input_vec2d);
      

    input_general(input_vec2d);

}
  
int input_nlopt(vector<vector<string> >& input_vec2d){
  
  
  

      
    get_vector_double(myinput.npar,"NLOPT_xl",input_vec2d,myinput.NLOPT_xl);
    get_vector_double(myinput.npar,"NLOPT_xu",input_vec2d,myinput.NLOPT_xu);
    get_vector_double(myinput.npar,"NLOPT_x0",input_vec2d,myinput.NLOPT_x0);
    getkey_int("NLOPT_nval",input_vec2d,myinput.NLOPT_nval);
    getkey_str("NLOPT_ALG",input_vec2d,myinput.NLOPT_ALG);
    getkey_double("NLOPT_etol", input_vec2d,myinput.NLOPT_etol);
}  

int input_geom_bond(vector<vector<string> >& input_vec2d){
  
  
  
  
getkey_int("geom_ngroup", input_vec2d,myinput.geom_ngroup);
getkey_int("geom_n", input_vec2d,myinput.geom_n);
vector<vector<double> > atom_vec2d;
getatom(myinput.geom_n,input_vec2d, atom_vec2d);
for(int i=0;i<myinput.geom_n;i++){
    for(int j=0;j<3;j++)myinput.geom_atom[i][j]=atom_vec2d[i][j];
}
vector<char> geom_group;
get_geom_group(myinput.geom_n ,input_vec2d ,geom_group);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_group[i]=geom_group[i];

vector<int> geom_xi;
get_geom_xi(myinput.geom_n,input_vec2d,geom_xi);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_xi[i]=geom_xi[i];
}  

int input_geom_3d(vector<vector<string> >& input_vec2d){

  
  
  
  
for(int i=0;i<50;i++)cout<<"#";
cout<<"In input_geom_3d"<<endl;
getkey_int("geom_ngroup", input_vec2d,myinput.geom_ngroup);
getkey_int("geom_n", input_vec2d,myinput.geom_n);
vector<vector<double> > atom_vec2d;
getatom(myinput.geom_n,input_vec2d, atom_vec2d);
for(int i=0;i<myinput.geom_n;i++){
    for(int j=0;j<3;j++)myinput.geom_atom[i][j]=atom_vec2d[i][j];
}
vector<char> geom_group;
get_geom_group(myinput.geom_n ,input_vec2d ,geom_group);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_group[i]=geom_group[i];
  
vector<int> geom_r;
get_geom_xngro("geom_r",myinput.geom_n,input_vec2d,geom_r);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_r[i]=geom_r[i];
vector<int> geom_theta;
get_geom_xngro("geom_theta",myinput.geom_n,input_vec2d,geom_theta);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_theta[i]=geom_theta[i];
vector<int> geom_phi;
get_geom_xngro("geom_phi",myinput.geom_n,input_vec2d,geom_phi);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_phi[i]=geom_phi[i];

  

}  



  
int input_general(vector<vector<string> >& input_vec2d){
  
  
  
  
  
for(int i=0;i<50;i++)cout<<"#";
cout<<"In input_general"<<endl;
int npot=0;
getkey_int("npot", input_vec2d,npot);
getpot(npot,input_vec2d,myinput.potmap);
  
getkey_int("ip_num", input_vec2d,myinput.ip_num );  
myinput.nstr_ALG="COBYLA";  
getkey_str("nstr_ALG",input_vec2d,myinput.nstr_ALG);  
getkey_int("nstr_num", input_vec2d,myinput.nstr_num);
  
#ifdef MINGW
memcpy(myinput.feffbat,"feffwin.Bat",CHAR_LEN);
memcpy(myinput.feffrun,"feffwin.Bat 1>1.log 2>2.log",CHAR_LEN);
memcpy(myinput.fdmneswin,"fdmnes_win32.exe 1>1.log 2>2.log",CHAR_LEN);
#else
memcpy(myinput.feffbat,"feff.sh",CHAR_LEN);
memcpy(myinput.feffrun,"feff.sh >& 1.log",CHAR_LEN);
memcpy(myinput.fdmneswin,"fdmnes >& 1.log ",CHAR_LEN);
#endif   
}

  
int input_feffprm(vector<vector<string> >& input_vec2d){
  
  
  
  
cout<<"In input_feffprm"<<endl;
    string strin;
    char*  cstrin;
    double ini,low,up;
    double din;
    int    iin;
  
    getkey_str("indir",input_vec2d,strin);
    strcpy(myinput.indir,strin.c_str());
  
    getkey_str("pdbname",input_vec2d,strin);
    cstrin=new char[CHAR_LEN];  
    strcpy(cstrin,strin.c_str());
  
  
  
    constr(myinput.indir,cstrin,myinput.pdbname);
  
   
  
    getkey_str("prmname",input_vec2d,strin);
    strcpy(cstrin,strin.c_str());
    constr(myinput.indir,cstrin,myinput.prmname);
      
  
    getkey_double3("es",input_vec2d,myinput.ini_es,myinput.lo_es,myinput.up_es);
    getkey_double3("ra",input_vec2d,myinput.ini_ra,myinput.lo_ra,myinput.up_ra);
}  

int input_ff2x(vector<vector<string> >& input_vec2d){
      
    if(myinput.runtype==5)
        {
        getkey_double3("ef_feff",input_vec2d,myinput.ini_ef_feff,myinput.lo_ef_feff,myinput.up_ef_feff);
        getkey_double3("conv_feff",input_vec2d,myinput.ini_conv_feff,myinput.lo_conv_feff,myinput.up_conv_feff);
        }
    if(myinput.runtype==6)
        {
        getkey_double("ef_feff_prm", input_vec2d,myinput.ef_feff_prm );
        getkey_double("conv_feff_prm", input_vec2d,myinput.conv_feff_prm);
        }
}  

int input_feffdiff(vector<vector<string> >& input_vec2d){


  
  
  

      
    string item;
    string strin;
    char*  cstrin;
    cstrin=new char[CHAR_LEN];  
    double ini,low,up;
    double din;
    int    iin;
    getkey_str("indir",input_vec2d,strin);
    strcpy(myinput.indir,strin.c_str());
  
    getkey_str("diff_exp_char",input_vec2d,strin);
    strcpy(cstrin,strin.c_str());
    constr(myinput.indir,cstrin,myinput.diff_exp_char);

  
  
    getkey_str("pdbname",input_vec2d,strin);
      
    strcpy(cstrin,strin.c_str());
  
  
  
constr(myinput.indir,cstrin,myinput.pdbname);
  


   
    getkey_str("prm_the_char",input_vec2d,strin);
  
    strcpy(cstrin,strin.c_str());
    constr(myinput.indir,cstrin,myinput.prm_the_char);
  

    getkey_double3("ex_diff",input_vec2d,myinput.ini_ex_diff,myinput.lo_ex_diff,myinput.up_ex_diff);
    

    getkey_double3("es_diff",input_vec2d,myinput.ini_es_diff,myinput.lo_es_diff,myinput.up_es_diff);
     
  
    getkey_double("es_prm", input_vec2d,myinput.es_prm );
    getkey_double("ra_prm", input_vec2d,myinput.ra_prm );
  

      

}


int input_fdmprm(vector<vector<string> >&input_vec2d){
    cout<<"In input_fdmprm"<<endl;
    string strin;
    char*  cstrin;
    double ini,low,up;
    double din;
    int    iin;

    getkey_int("nopt_conv", input_vec2d,myinput.nopt_conv);  
    myinput.n_conv=5;
    getkey_int("n_conv", input_vec2d,myinput.n_conv);  

    get_vector_double(myinput.n_conv,"fdm_up_conv",input_vec2d,myinput.fdm_up_conv);
    get_vector_double(myinput.n_conv,"fdm_low_conv",input_vec2d,myinput.fdm_low_conv);
    get_vector_double(myinput.n_conv,"fdm_x0_conv",input_vec2d,myinput.fdm_x0_conv);
}

  




  
  
  
  
int getkey_str(const string a,  vector<vector<string> >& input_vec2d,string& key_str){
    int find=0;  
      
if(myinput.display==1) cout<<"we are searing "<<a<<endl;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==a){
                 
                string ss=*(j+1);
                key_str=ss;
                 
                find=1;
                if(myinput.display==1) cout<<a<<" is "<<key_str<<endl;
        }  
        }  
         
  
  
  
  

    return 0;
}  


int getkey_double(const string a,  vector<vector<string> >& input_vec2d,double& key_double){
    if(myinput.display==1)cout<<"We are searing "<<a<<endl;
      
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==a){
                string ss=*(j+1);
                key_double=atof(ss.c_str());
                if(myinput.display==1) cout<<a<<" is "<<key_double<<endl;
        }  
        }  
    return 0;
}  
int getkey_int(const string a,  vector<vector<string> >& input_vec2d,int& key_int){
    if(a!="display"&&myinput.display==1) cout<<"We are searching "<<a<<endl;
      
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==a){
                string ss=*(j+1);

                key_int=atoi(ss.c_str());
                if(a!="display"&&myinput.display==1)cout<<a<<" is "<<key_int<<endl;
        }  
        }  
    return 0;
}  
  
int getkey_double3(const string a,  vector<vector<string> >& input_vec2d ,double& key_1,double & key_2,double & key_3){
        if(myinput.display==1) cout<<"We are searching "<<a<<endl;
          
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==a){
                string ss1=*(j+1);
                key_1=atof(ss1.c_str());
                string ss2=*(j+2);
                key_2=atof(ss2.c_str());
                string ss3=*(j+3);
                key_3=atof(ss3.c_str());
                if(myinput.display==1)cout<<a<<" ini low_bind up_bond are "<<key_1<<" "<<key_2<<" "<<key_3<<endl;
        }  
        }  
    return 0;
}  
  
int get_vector_double(int n,const string vecname, vector<vector<string> >& input_vec2d,vector<double>& vec_out){
      
          
        if(myinput.display==1)cout<<"We are searching "<<vecname<<endl;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==vecname)iblock=i;
        }  
        vector<string> item=*iblock;
        if(myinput.display==1)cout<<vecname<<" are ";
        for(int ip=1;ip<n+1;++ip){

            string val_str=item[ip];
            double val=atof(val_str.c_str());
            vec_out.push_back(val);
            if(myinput.display==1) cout<<" "<<val;
        }
        if(myinput.display==1)cout<<endl;

    return 0;
}  
  
int getpot(int npot, vector<vector<string> >& input_vec2d,map<string,int>& potmap){
    if(myinput.display==1) cout<<"We are searing pot"<<endl;
      
          
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j=="pot")iblock=i+1;
        }  
        for(int ip=0;ip<npot;++ip){
            vector<string> item=*(iblock+ip);
            string ele=item[0];
            int pot=atoi(item[1].c_str());
            potmap[ele]=pot;
            if(myinput.display==1)cout<<ele<<" "<<potmap[ele]<<endl;
        }

    return 0;
}  
  
int getatom(int nchange, vector<vector<string> >& input_vec2d,vector<vector<double> >& atom_vec2d){
      
         
         
       if(myinput.display==1) cout<<"We are searching representative atom "<<endl;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();

        if(*j=="atom"){iblock=i+1;}
          
        }  
          
        for(int ip=0;ip<nchange;++ip){
            vector<string> item=*(iblock+ip);


            double x=atof(item[0].c_str());
            double y=atof(item[1].c_str());
            double z=atof(item[2].c_str());
            vector<double> atom_vec;
            atom_vec.push_back(x);
            atom_vec.push_back(y);
            atom_vec.push_back(z);
            atom_vec2d.push_back(atom_vec);
            if(myinput.display==1)cout<<x<<" "<<y<<" "<<z<<endl;
        }

    return 0;
}  

int get_geom_group(int geom_n, vector<vector<string> >& input_vec2d,vector<char>& geom_group){
    if(myinput.display==1) cout<<"We are searching atom group char"<<endl;
      
          
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j=="geom_group")iblock=i;
        }  
        vector<string> item=*iblock;
        for(int ip=1;ip<geom_n+1;++ip){

            char group_char=item[ip][0];
            geom_group.push_back(group_char);
            if(myinput.display==1) cout<<" "<<group_char;
        }
        if(myinput.display==1) cout<<endl;

    return 0;
}  

int get_geom_xi(int geom_n, vector<vector<string> >& input_vec2d,vector<int>& geom_xi){
    if(myinput.display==1)cout<<"We are searching xi of group"<<endl;
      
          
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j=="geom_xi")iblock=i;
        }  
        vector<string> item=*iblock;
        for(int ip=1;ip<geom_n+1;++ip){

            int xi=atoi(item[ip].c_str());
            geom_xi.push_back(xi);
            if(myinput.display==1)cout<<" "<<xi;
        }
        if(myinput.display==1)cout<<endl;

    return 0;
}  
  
int get_geom_xngro(const string key_char,int geom_n, vector<vector<string> >& input_vec2d,vector<int>& geom_xi){
    if(myinput.display==1)cout<<"We are searching xi of "<<key_char<<endl;
      
          
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==key_char)iblock=i;
        }  
        vector<string> item=*iblock;
        for(int ip=1;ip<geom_n+1;++ip){

            int xi=atoi(item[ip].c_str());
            geom_xi.push_back(xi);
            if(myinput.display==1)cout<<" "<<xi;
        }
        if(myinput.display==1)cout<<endl;

    return 0;
}  











  
int readinp(const char* finp,vector<vector<string> >& input_vec2d){
    ifstream fin(finp);
    string s;
    while( getline(fin,s) )  
    {   vector<string> oneline;
        if(s.length()==0) continue;  
         
        boost::char_separator<char> mysep("  ");  
        boost::tokenizer< boost::char_separator<char> > mytok(s,mysep);  
        int i_tok=0;
        vector<string> mytok_vec;

          
        for(boost::tokenizer<boost::char_separator<char> >::iterator i=mytok.begin();i!=mytok.end();++i)
        {i_tok=i_tok+1;
        mytok_vec.push_back(*i);
        };
          
        if(mytok_vec[0][0]=='#')continue;  
        input_vec2d.push_back(mytok_vec);
    }  
return 0;
}

int dis_vec2d(vector<vector<string> > input_vec2d){
    int idx=0;
for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
    vector<string> vv=*i;
    std::cout<<"line "<<idx<<" ";
    for(vector<string>::iterator j=vv.begin();j!=vv.end();++j){
            std::cout<<*j<<" ";

    }  
    std::cout<<std::endl;
    ++idx;
}  
return 0;
}  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

#endif   
