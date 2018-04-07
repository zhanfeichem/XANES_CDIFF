#ifndef DL_FDM
#define DL_FDM
using namespace std;
//CPP header
#include <iostream>
#include <fstream>
#include<sstream>
#include<vector>
#include<map>
#include <string>
//boost header
#include<boost/tokenizer.hpp>//分割字符串
//#include<boost/lexical_cast.hpp>//转换类型
#include<boost/format.hpp>//格式化输出
//C header
#include<stdlib.h>//使用exit()
#include<math.h>//使用三角函数
#include <string.h>//c风格字符串操作
#include<time.h> //时间做随机数种子
#include<stdlib.h>
#include <unistd.h>//切换工作目录 _acess函数
#include <setjmp.h>
#include<signal.h>
#include <sys/stat.h>//mkdir函数
#include <sys/types.h>//mkdir函数
//#include<sys/wait.h>
//My header


#include "subfdmdiff.hpp"//差分普 暂时屏蔽
using namespace fdmdiff;
#include "subfdm.hpp"
using namespace fdm;
//
#include "mytool.hpp"
#include "glo.hpp"
#include "read_dat.hpp"//read_fdm
double dl_mycfit_fdm(int ,double* );
//int fdminp(vector<PDB> pdbline);
//int fdmnesinp(vector<PDB> pdbline,char workdir_char[CHAR_LEN]);//使用runfdm中的版本
double dl_save(int npar,double* x);

double dl_mycfit_fdm(int npar,double* x){
//char mycom[CHAR_LEN];//所有shell命令用此执行 记得每次清空
double result = 1e21;
vector<PDB> pdbline;
getpdb(pdbline);

//改变几何参数
change_shell_ini(pdbline);

//int n_change=1;
//change_shell(n_change,'B',x[0],pdbline);
//change_shell(n_change,'C',x[1],pdbline);
//change_shell(n_change,'D',x[1],pdbline);
//change_shell(n_change,'E',x[2],pdbline);
//change_shell(n_change,'F',x[2],pdbline);
//改变几何参数

//double atom2[3]={0.859,-1.195,1.404};
//double atom3[3]={0.275,-1.726,-1.219};
//double atom4[3]={-1.949,-0.467,0.667};
//double atom5[3]={1.918,0.774,-0.430};
//double atom6[3]={-0.053,1.489,1.507};
//int n_change=0;
//double chang1[3]={x[0],0.0,0.0};
//double chang2[3]={x[1],0.0,0.0};
//double chang3[3]={x[2],0.0,0.0};
////bpy axis
//change_lig(n_change,'B',chang1,atom2,pdbline);
////bpy plane
//change_lig(n_change,'C',chang2,atom3,pdbline);
////dpa
//change_lig(n_change,'D',chang3,atom4,pdbline);
//change_lig(n_change,'E',chang3,atom5,pdbline);
//change_lig(n_change,'F',chang3,atom6,pdbline);
if(myinput.geom_model==1){
int n_change=myinput.geom_ngroup;

for(int i=0;i<myinput.geom_n;i++){
        int nn=myinput.geom_xi[i];
        double change[3]={x[nn],0,0};
        double atom[3]={myinput.geom_atom[i][0],myinput.geom_atom[i][1],myinput.geom_atom[i][2]};
        change_lig(n_change,myinput.geom_group[i],change,atom,pdbline);
}//end for

}//end if
else if(myinput.geom_model==11){
int n_change=myinput.geom_ngroup;

for(int i=0;i<myinput.geom_n;i++){

        int nn=myinput.geom_xi[i];

        double change[3]={x[nn],0,0};

        //cout<<"dr  "<<change[0]<<" dtheta "<<change[1]<<" dphi "<<change[2]<<endl;
        double atom[3]={myinput.geom_atom[i][0],myinput.geom_atom[i][1],myinput.geom_atom[i][2]};
        change_lig_val(n_change,myinput.geom_group[i],change,atom,pdbline);
}//end for

}//end else11

else if(myinput.geom_model==2){
int n_change=myinput.geom_ngroup;
for(int i=0;i<myinput.geom_n;i++){
        int ri=myinput.geom_r[i];
        int thetai=myinput.geom_theta[i];
        int phii=myinput.geom_phi[i];
        //
        double change[3]={x[ri],x[thetai],x[phii]};
        //cout<<"r theta phi"<<x[ri]<<x[thetai]<<x[phii]<<endl;
        //cout<<"dr  "<<change[0]<<" dtheta "<<change[1]<<" dphi "<<change[2]<<endl;
        double atom[3]={myinput.geom_atom[i][0],myinput.geom_atom[i][1],myinput.geom_atom[i][2]};
        change_lig_val(n_change,myinput.geom_group[i],change,atom,pdbline);
}//end for

}//end else





//输出pdb
save_pdb(pdbline,npar,x);
//建立文件夹
char workdir_char[CHAR_LEN];
//myworkdir(workdir_char);//获得随机数字符串 作为运行目录
//fdmnes使用固定文件夹
strcpy(workdir_char,"fdmnesrun");
//feffinp
//feffinp(pdbline);
fdmnesinp(pdbline,workdir_char);//工作目录下建立feffinp
//debug

//运行feff
//兼容问题现在不更改目录 采取拷贝运行脚本到指定文件夹  遗留问题fdmnes如何处理？？
char dirnow_char[CHAR_LEN];//保存当前目录以便返回
getcwd(dirnow_char,sizeof(dirnow_char));
chdir(workdir_char);//如果改变目录在if setjnp之前则无法删除工作目录 因为使用的是相对路径 当前目录不对

//复制fdmnes需要文件  复制文件函数

if(myinput.runCore==1)system(myinput.fdmneswin);


//运行部分结束
chdir(dirnow_char);//从运行目录返回主目录
//判断xmu.dat是否存在
char xmupos_char[CHAR_LEN];
memcpy(xmupos_char,workdir_char,CHAR_LEN);
strcat(xmupos_char,"/out.txt");//输出文件名
if(access(xmupos_char,0)==-1) {

//        char rm_char[CHAR_LEN]=RMDIR;//given blank in the edn of the commad
//
//        strcat(rm_char,workdir_char);
//
//        system(rm_char);
        result=1e22;
        wrongfile("saveall",x,npar,"noxmu");
        remove("./fdmnesrun/out.bk");
        remove("./fdmnesrun/out_conv.bk");
        return result;
}//if block ends

//几何结构变量传到全局变量mypar结构体
myinput.fdm_ngeom=npar;
//myinput.pargeom=(double *)malloc(mypar.ngeom*sizeof(double));
myinput.fdm_pargeom=vector<double>(npar);
for(int i=0;i<npar;i++)myinput.fdm_pargeom[i]=x[i];

//保存训练数据集
dl_save(npar,x);
//
remove("./fdmnesrun/out.bk");
remove("./fdmnesrun/out_conv.bk");
rename("./fdmnesrun/out.txt","./fdmnesrun/out.bk");
rename("./fdmnesrun/out_conv.txt","./fdmnesrun/out_conv.bk");
//return result;
return 0;


}//end of mycfit函数结束




double dl_save(int npar,double* x){
    vector<vector<double> > dat_xmu;
    vector<double> x_xmu,y_xmu;
    //从文件读取数据
    read_dat_fdm("fdmnesrun/out_conv.txt",2,dat_xmu);//打开卷积谱
    x_xmu=dat_xmu[0];y_xmu=dat_xmu[1];

    //数据保存
    ofstream dlsave("dl_spc.txt",ios::app);
    for(unsigned i=0;i<=y_xmu.size()-1;++i)
       dlsave<<boost::format("%15.8f")%(y_xmu[i]);
    dlsave<<endl;//此条谱结束
    //输出参数
    ofstream dlpar("dl_par.txt",ios::app);
    for(unsigned i=0;i<npar;i++){
        dlpar<<x[i]<<" ";
        std::cout<<x[i]<<"  ";
    }
    dlpar<<endl;
    cout<<endl;

}













#endif // DL_FDM
