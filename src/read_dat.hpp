#ifndef READ_DAT
#define READ_DAT


#include <iostream>
#include<fstream>
#include<vector>
#include<boost/tokenizer.hpp> 
 
#include<stdlib.h>
#include<boost/format.hpp> 
using namespace std;
int read_dat(const char* fn,int n_col,vector<vector<double> > & vec_out)
{    
     

    ifstream mu(fn);
    string oneline;
    vector<string> vec_xmu; 
    vector<vector<string> > vec_tot; 
    string xmu_comment("#"); 
     
    while(getline(mu,oneline)){

          string ss;
          ss.push_back(oneline[0]);
          if (ss.compare(string("#"))==0)continue;
          vec_xmu.push_back(oneline);
          } 


     
    for(vector<string>::iterator i=vec_xmu.begin();i!=vec_xmu.end();++i){
        boost::char_separator<char> mysep("  "); 
        boost::tokenizer< boost::char_separator<char> > mytok(*i,mysep);
        vector<string> mytok_vec;
        for(boost::tokenizer<boost::char_separator<char> >::iterator i=mytok.begin();i!=mytok.end();++i){
            mytok_vec.push_back(*i);

        } 

        vec_tot.push_back(mytok_vec);

    } 



     
    vector<vector<double> > dat; 
    for(int i=0;i<=n_col-1;++i){
        vector<double> vv;
        dat.push_back(vv);
    } 

     
    for(vector<vector<string> >::iterator i=vec_tot.begin();i!=vec_tot.end();++i){
        int n_dat =0;
        for(vector<string>::iterator j=i->begin();j!=i->end();++j){
            double dd = 0;
             
            dd=atof((*j).c_str());
            dat[n_dat].push_back(dd);
            ++n_dat;
             
        } 

    } 
vec_out=dat;
return 0;} 


int read_dat_fdm(const char* fn,int n_col,vector<vector<double> > & vec_out)
{    
     
 
 
 
 



    ifstream mu(fn);
    string oneline;
    vector<string> vec_xmu; 
    vector<vector<string> > vec_tot; 
    string xmu_comment("#"); 
     
    int flag=0;
    while(getline(mu,oneline)){


          string ss;
          ss.push_back(oneline[0]);
          if (ss.compare(string("#"))==0)continue;

          if (oneline.find("xanes",0)!=string::npos){flag=1;continue;} 

          if (flag==0) continue;
          vec_xmu.push_back(oneline);

    
          } 
    


     
    for(vector<string>::iterator i=vec_xmu.begin();i!=vec_xmu.end();++i){
        boost::char_separator<char> mysep("  "); 
        boost::tokenizer< boost::char_separator<char> > mytok(*i,mysep);
        vector<string> mytok_vec;
        for(boost::tokenizer<boost::char_separator<char> >::iterator i=mytok.begin();i!=mytok.end();++i){
            mytok_vec.push_back(*i);

        } 

        vec_tot.push_back(mytok_vec);

    } 



     
    vector<vector<double> > dat; 
    for(int i=0;i<=n_col-1;++i){
        vector<double> vv;
        dat.push_back(vv);
    } 

     
    for(vector<vector<string> >::iterator i=vec_tot.begin();i!=vec_tot.end();++i){
        int n_dat =0;
        for(vector<string>::iterator j=i->begin();j!=i->end();++j){
            double dd = 0;
             
            dd=atof((*j).c_str());
            dat[n_dat].push_back(dd);
            ++n_dat;
             
        } 

    } 
vec_out=dat;
return 0;} 
#endif
