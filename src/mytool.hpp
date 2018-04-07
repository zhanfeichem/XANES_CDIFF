#ifndef MYTOOL
#define MYTOOL










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
 

 
#include "glo.hpp"

struct PDB{
int pot;
string comment;
string ele;
string group;
double x,y,z;
double r,theta,phi;
double xn,yn,zn;
double rn,thetan,phin;
};

void array2string(double* x,int nx,string& strout);
void wrongfile(char* folder,double* x,int nx,char* errorchar);



int getpdb(vector<PDB>& pdbline_out);
int save_pdb(vector<PDB> pdbline,int npar,double* x);
int change_shell_ini(vector<PDB>& pdbline);
int change_shell(int n_group,char group_char,double rt,vector<PDB>& pdbline);
int change_lig(int n_group,char group_char,double change_percent[3],double atom[3],vector<PDB>& pdbline);
int myworkdir(char dir_char[CHAR_LEN]);
int car2sph(double x,double y,double z,vector<double>& re);
int sph2car(double r,double theta,double phi,vector<double>& re);
int my_itoa(int val, char* buf);

void wrongfile(char* folder,double* x,int nx,char* errorchar){
     
    string xstring;
    array2string(x,nx,xstring);
    ostringstream ostr;
    ostr<<folder<<"/";
    ostr<<xstring;
    ostr<<errorchar<<".txt";
    string ostr_str=ostr.str();


    ofstream mysave(ostr_str.c_str());
}
void array2string(double* x,int nx,string& strout){
 
    ostringstream ostr;
    for(int i=0;i<nx;++i)ostr<<x[i]<<"_";
    strout=ostr.str();
}

int getpdb(vector<PDB>& pdbline_out){
     
    ifstream fin(myinput.pdbname);
    map<string,int> potmap=myinput.potmap;
 
 
 
 
 
     
    vector<string> myline; 
    vector<PDB>    pdbline; 
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
        if(i_tok!=8){
            cerr<<"number of PDB item is wrong "<<i_tok<<endl;
            cerr<<"item1  "<<mytok_vec[0]<<"  item2  "<<mytok_vec[1]<<"  item3  "<<mytok_vec[2]<<"  item4  "<<mytok_vec[3]<<"  item5  "<<mytok_vec[4]<<"  item6 x "<<mytok_vec[5]<<"  item7 y  "<<mytok_vec[6]<<"  item8z  "<<mytok_vec[7]<<endl;
            exit(1);
        };

       PDB ipdb;
 
 
 
       ipdb.x=atof(mytok_vec[5].c_str());
       ipdb.y=atof(mytok_vec[6].c_str());
       ipdb.z=atof(mytok_vec[7].c_str());
       ipdb.group=mytok_vec[3];
       ipdb.ele=mytok_vec[2];
       pdbline.push_back(ipdb);
       myline.push_back(s);
    }; 
     
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
        vector<double> re;
        car2sph((i->x),(i->y),(i->z),re);
        (i->r)=re[0];
        (i->theta)=re[1];
        (i->phi)=re[2];
    } 
     
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
        if(potmap.count((i->ele))==0){
            cout<<(i->ele)<<"do not in potmap"<<endl;exit(-1);
        };
        (i->pot)=potmap[(i->ele)];
    } 
    pdbline_out=pdbline; 
return 0;
}; 
int save_pdb(vector<PDB> pdbline,int npar,double* x){
    ostringstream oo;
    oo<<"savepdb/";
    for(int i=0;i<=npar-1;++i){
        oo<<x[i]<<"_";
    }
    oo<<"savepdb.pdb";
    string ss;ss=oo.str();
    const char* cc=ss.c_str();
    ofstream outpdb(cc); 
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i)
        outpdb<< boost::format("ATOM%7.0f%3s%6s%6.0f%4s%8.3f%8.3f%8.3f\n")%1%(i->ele)%(i->group)%1%(" ")%(i->xn)%(i->yn)%(i->zn)<<endl;
return 0;} 
 
int change_shell_ini(vector<PDB>& pdbline){
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
        (i->rn)=(i->r);
        (i->thetan)=(i->theta);
        (i->phin)=(i->phi);
    } 
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
        vector<double> re;
        sph2car(i->rn,i->thetan,i->phin,re);
        (i->xn)=re[0];
        (i->yn)=re[1];
        (i->zn)=re[2];
    }
return 0;} 
int change_shell(int n_group,char group_char,double rt,vector<PDB>& pdbline){
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
            if(i->group[n_group]==group_char){

                (i->rn)=(i->r)+rt*(i->r);
                (i->thetan)=(i->theta);
                (i->phin)=(i->phi);
            } 

    } 
     
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
        if(i->group[n_group]==group_char){
            vector<double> re;
            sph2car(i->rn,i->thetan,i->phin,re);
            (i->xn)=re[0];
            (i->yn)=re[1];
            (i->zn)=re[2];
        } 
   } 
return 0;} 



 
int change_lig(int n_group,char group_char,double change_percent[3],double atom[3],vector<PDB>& pdbline){
     
     
    vector<double> atom_sph;
    vector<double> delta; 
    car2sph(atom[0],atom[1],atom[2],atom_sph);
    for(int i=0;i<3;i++)atom_sph[i]=atom_sph[i]+change_percent[i]*atom_sph[i];
    sph2car(atom_sph[0],atom_sph[1],atom_sph[2],delta);
    for(int i=0;i<3;i++)delta[i]=delta[i]-atom[i];




    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){

            if(i->group[n_group]==group_char){

                (i->xn)=(i->x)+delta[0];
                (i->yn)=(i->y)+delta[1];
                (i->zn)=(i->z)+delta[2];
            } 

    } 
     
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
        if(i->group[n_group]==group_char){
            vector<double> re;
            car2sph(i->xn,i->yn,i->zn,re);
            (i->rn)=re[0];
            (i->thetan)=re[1];
            (i->phin)=re[2];
        } 
   } 
return 0;} 

 
int change_lig_val(int n_group,char group_char,double change_val[3],double atom[3],vector<PDB>& pdbline){
     
     
    vector<double> atom_sph;
    vector<double> delta; 
    car2sph(atom[0],atom[1],atom[2],atom_sph);
     
    for(int i=0;i<3;i++)atom_sph[i]=atom_sph[i]+change_val[i]; 
    sph2car(atom_sph[0],atom_sph[1],atom_sph[2],delta);
    for(int i=0;i<3;i++)delta[i]=delta[i]-atom[i];




    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){

            if(i->group[n_group]==group_char){

                (i->xn)=(i->x)+delta[0];
                (i->yn)=(i->y)+delta[1];
                (i->zn)=(i->z)+delta[2];
            } 

    } 
     
    for(vector<PDB>::iterator i=pdbline.begin();i!=pdbline.end();++i){
        if(i->group[n_group]==group_char){
            vector<double> re;
            car2sph(i->xn,i->yn,i->zn,re);
            (i->rn)=re[0];
            (i->thetan)=re[1];
            (i->phin)=re[2];
        } 
   } 
return 0;} 




int myworkdir(char dir_char[CHAR_LEN]){
     
     

    srand(time(0));
    int ii=rand();
    char ss[CHAR_LEN];
     
    my_itoa(ii,ss);
     
   while(access(ss,0) != -1){
       ii=rand();
        
       my_itoa(ii,ss);
   } 
    #ifdef MINGW
    mkdir(ss);
    #else
    mkdir(ss,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #endif  
    memcpy(dir_char,ss,CHAR_LEN); 

return 0;} 
int car2sph(double x,double y,double z,vector<double>& re){
re.clear();
double r,theta,phi = 0;
double z2 = 0;
z2=pow(x,2)+pow(y,2)+pow(z,2);
r=sqrt(z2);
theta=acos(z/r);
phi=atan2(y,x);
re.push_back(r);
re.push_back(theta);
re.push_back(phi);
 
return 0;
};
int sph2car(double r,double theta,double phi,vector<double>& re){
re.clear();
if(r==0){re.push_back(0);re.push_back(0);re.push_back(0);}
double x,y,z = 0;
x=r*sin(theta)*cos(phi);
y=r*sin(theta)*sin(phi);
z=r*cos(theta);
re.push_back(x);
re.push_back(y);
re.push_back(z);
 
return 0;
};
int my_itoa(int val, char* buf){
    const unsigned int radix = 10;

    char* p;
    unsigned int a;         
    int len;
    char* b;             
    char temp;
    unsigned int u;

    p = buf;

    if (val < 0)
    {
        *p++ = '-';
        val = 0 - val;
    }
    u = (unsigned int)val;

    b = p;

    do
    {
        a = u % radix;
        u /= radix;

        *p++ = a + '0';

    } while (u > 0);

    len = (int)(p - buf);

    *p-- = 0;

     
    do
    {
        temp = *p;
        *p = *b;
        *b = temp;
        --p;
        ++b;

    } while (b < p);

    return len;
}

int print_atom(double* atom){
    cout<<atom[0]<<" "<<atom[1]<<" "<<atom[2]<<endl;
return 0;
} 


#endif  
