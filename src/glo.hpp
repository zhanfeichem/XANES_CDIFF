#ifndef GLO
#define GLO

//处理兼容性
//#define CENTOS
#define MINGW
//#define USENOMAD//由于有的系统不兼容NOMAD 只有定义了USENOMAD才包含NOMAD模块内容
//兼容性选型结束

#ifdef CYGWIN
#include<io.h>
#endif // CYGWIN
#ifdef CENTOS
#include <sys/io.h>//_access函数判断文件存在
#endif // CENTOS

//shell命令宏定义 cygwin 和 linux一致
#ifdef MINGW//注意在mingw32下 wait和signal不能使用 mkdir的参数只有一个和cygwin不同
#define RMDIR "rd /s/q "//注意命令后都空格
//#include<sys/wait.h>//没有wait.h wait和signal不能使用
#include<io.h>
#else
#define RMDIR "rm -rf "//pay attention there should be a space after command
//#define CHMOD "chmod 777 "
//#define CP "cp "
#endif // MINGW
//C语言库
#include<string.h>
//STL
#include<map>
#include<vector>
#include<fstream>
using namespace std;
//boost
#include<boost/tokenizer.hpp>//分割字符串
//#include<boost/lexical_cast.hpp>//转换类型
#include<boost/format.hpp>//格式化输出
//
#include<cstdlib>

#define TIME_OUT 1200// limit of time exceed which can cause the alarm functon
#define CHAR_LEN 50//len of char array
//#define DIS// 定义是否开启显示 因为nomad的使用禁止程序显示

//NLOPT_LN_COBYLA
//NLOPT_LN_BOBYQA
//NLOPT_LN_NEWUOA_BOUND
//NLOPT_LN_PRAXIS
//NLOPT_LN_NELDERMEAD
//NLOPT_LN_SBPLX


#define NLOPT_CONV NLOPT_LN_BOBYQA
#define NLOPT_SUB  NLOPT_LN_COBYLA
//调试宏选项
//#define RUNOFF//关闭feff脚本的运行
//int print_atom1(double* atom);
//几何结构变化的相关宏
#define GEOM_N 20
//
struct par_str{
    int ngeom;
    double* pargeom;
    int nconv;
    double* parconv;
    double* parsub;
};//传递各层参数 geom结构参数 conv卷积参数 sub ratio和energyshift参数
struct Myin{
    //控制类参数
    int display;
    int runtype;
    int npar;
    int runCore;//feffon 0 无计算 1 prm 2 diff
    int savedir;
    //设置目标函数
    int OBJ;//1 默认 2 R因子 weight=1
    double OBJ_eps;//设置R因子中的eps为常数
    double OBJ_afit,OBJ_bfit;//设置目标函数使用的实验谱能量区间
    double OBJ_aout,OBJ_bout;//设置输出区间
    //NLOPT 参数
    string OPT;
    string NLOPT_ALG;
    int NLOPT_nval;
    double NLOPT_etol;
    vector<double> NLOPT_x0;
    vector<double> NLOPT_xl;
    vector<double> NLOPT_xu;
    //控制类参数结束
    char indir[CHAR_LEN];//输入文件夹

    char feffbat[CHAR_LEN];
    char feffrun[CHAR_LEN];
    char fdmneswin[CHAR_LEN];
    char pdbname[CHAR_LEN];
    char prmname[CHAR_LEN];//基态实验谱的名称
    //差分拟合相关文件
    char prm_the_char[CHAR_LEN];//基态理论谱的名称
    char diff_exp_char[CHAR_LEN];//

    int  ip_num;//插值点数量
    string nstr_ALG;//内层优化算法
    int  nstr_num;//es和ra等非结构参数拟合次数

    double up_es,lo_es;//enery shift 和 ratio的上下限制和初始值
    double up_ra,lo_ra;
    double ini_es,ini_ra;
    double ini_ef_feff,lo_ef_feff,up_ef_feff;//feff efermi
    double ini_conv_feff,lo_conv_feff,up_conv_feff;//feff 展宽
    //差分谱拟合参数
    double es_prm;//已经拟合好的基态能移
    double ra_prm;//已经拟合好的归一化因子
    double ef_feff_prm,conv_feff_prm;//ff2x使用的基态费米面移动和展宽
    double up_es_diff,lo_es_diff;//相对能移
    double up_ex_diff,lo_ex_diff;//激发比例
    double ini_es_diff,ini_ex_diff;
    double** atom_array;
    map<string,int> potmap;
    //卷积拟合相关参数
    int nopt_conv;//卷积拟合的次数
    int n_conv;//一般为Efermi Gamma_hole Gamma_max Ecent Elarg 5个参数

//    double* up_conv;
//    double* low_conv;
//    double* x0_conv;
    vector<double> fdm_up_conv;
    vector<double> fdm_low_conv;
    vector<double> fdm_x0_conv;
    //fdmnes中间传递参数
    vector<double> fdm_parconv;
    int fdm_ngeom;
    vector<double> fdm_pargeom;


//参数传递结构体
par_str mypar;

//几何结构变化相关参数
int geom_model;//选择几何变化的类型
int geom_ngroup;//使用的几何变换分组字段第一组为0 类推
int geom_n;//几何变换次数
double geom_atom[GEOM_N][3];//代表原子坐标数组
char  geom_group[GEOM_N];
//键长模式使用
int  geom_xi[GEOM_N];
//键长角度模式
int geom_r[GEOM_N];
int geom_theta[GEOM_N];
int geom_phi[GEOM_N];
};//定义保存输入参数的结构体
//定义初始变量 之后转换为输入文件
//定义全局变量
Myin myinput;
par_str mypar;
int constr(char* s1,char* s2,char* out){
    //连接字符串s1+s2
    memcpy(out,s1,CHAR_LEN);
    strcat(out,s2);
return 0;
}//函数结束

int readinp(const char* ,vector<vector<string> >& );
int dis_vec2d(vector<vector<string> > );//显示两层字符串
//get keywords
int getkey_str(const string , vector<vector<string> >& ,string& );
int getkey_double(const string , vector<vector<string> >& ,double& );
int getkey_int(const string,  vector<vector<string> >& ,int& );
//
int getkey_double3(const string ,  vector<vector<string> >& ,double& ,double & ,double & );
int get_vector_double(int ,const string , vector<vector<string> >& ,vector<double>& );
//
int getpot(int, vector<vector<string> >& ,map<string,int>& );
//
int getatom(int,vector<vector<string> >& ,vector<vector<double> >& );
int get_geom_group(int , vector<vector<string> >& ,vector<char>& );
int get_geom_xi(int , vector<vector<string> >&,vector<int>& );
int get_geom_xngro(const string ,int , vector<vector<string> >& ,vector<int>& );
//input function
int input_general(vector<vector<string> >&);//通用输入
int input_nlopt(vector<vector<string> >&);
int input_feffdiff(vector<vector<string> >&);
int input_feffprm(vector<vector<string> >&);
int input_ff2x(vector<vector<string> >&);
int input_fdmprm(vector<vector<string> >&);


int input_geom_bond(vector<vector<string> >&);
int input_geom_3d(vector<vector<string> >&);
//
int inputfile(){
    string strin;

    const char* finp="input.inp";
    vector<vector<string> > input_vec2d;
    readinp(finp,input_vec2d);
    dis_vec2d(input_vec2d);//最开始显示读入的字符串

    //display变量第一个读入
    myinput.display=0;//设置初始不显示
    getkey_int("display",input_vec2d,myinput.display);
    //
    getkey_int("npar",input_vec2d,myinput.npar);
    myinput.runCore=1;//初始设置runCore
    getkey_int("runCore",input_vec2d,myinput.runCore);
    myinput.savedir=0;//初始设置savedir
    getkey_int("savedir",input_vec2d,myinput.savedir);


    //选择优化算法
    myinput.OPT="NOMAD";//默认设置NOMAD
    getkey_str("OPT",input_vec2d,myinput.OPT);
    if(myinput.OPT=="NLOPT"||myinput.OPT=="NOMADLIB") input_nlopt(input_vec2d);//读取NLOPT设置

    //目标函数设置
    //myinput.OBJ=1;//使用之前的目标函数
    getkey_int("OBJ",input_vec2d,myinput.OBJ);
    if(myinput.OBJ==2)getkey_double("OBJ_eps",input_vec2d,myinput.OBJ_eps);//获取R因子的epsilon
    if(myinput.OBJ==3)getkey_double("OBJ_afit",input_vec2d,myinput.OBJ_afit);//目标函数3 获取拟合区间ea
    if(myinput.OBJ==3)getkey_double("OBJ_bfit",input_vec2d,myinput.OBJ_bfit);//eb
    if(myinput.OBJ==3)getkey_double("OBJ_aout",input_vec2d,myinput.OBJ_aout);//输出区间 输出区间一定包含fit区间
    if(myinput.OBJ==3)getkey_double("OBJ_bout",input_vec2d,myinput.OBJ_bout);

    //

    getkey_int("runtype",input_vec2d,myinput.runtype);
    if(myinput.runtype==1||myinput.runtype==5) input_feffprm(input_vec2d);//5 展宽的FEFF拟合
    if(myinput.runtype==2||myinput.runtype==6) input_feffdiff(input_vec2d);//6 ff2x展宽的FEFF拟合
    if(myinput.runtype==5||myinput.runtype==6) input_ff2x(input_vec2d);//5 模式读取feff展宽参数
    //if(myinput.runtype==5) input_ff2x(input_vec2d);//6 模式读取feff展宽参数

    if(myinput.runtype==3) input_feffprm(input_vec2d);//先读取feff版的基态函数
    if(myinput.runtype==3) input_fdmprm(input_vec2d);//fdmnes的基态读取函数
    if(myinput.runtype==4) input_feffdiff(input_vec2d);
    if(myinput.runtype==4) input_fdmprm(input_vec2d);//FDM差分普也使用基态fdm的输入函数

    getkey_int("geom_model",input_vec2d,myinput.geom_model);
    if(myinput.geom_model==1||myinput.geom_model==11) input_geom_bond(input_vec2d);
    if(myinput.geom_model==2) input_geom_3d(input_vec2d);
    //通用读取部分

    input_general(input_vec2d);

}
//
int input_nlopt(vector<vector<string> >& input_vec2d){
//    const char* finp="input.inp";
//    vector<vector<string> > input_vec2d;
//    readinp(finp,input_vec2d);

    //read
    get_vector_double(myinput.npar,"NLOPT_xl",input_vec2d,myinput.NLOPT_xl);
    get_vector_double(myinput.npar,"NLOPT_xu",input_vec2d,myinput.NLOPT_xu);
    get_vector_double(myinput.npar,"NLOPT_x0",input_vec2d,myinput.NLOPT_x0);
    getkey_int("NLOPT_nval",input_vec2d,myinput.NLOPT_nval);
    getkey_str("NLOPT_ALG",input_vec2d,myinput.NLOPT_ALG);
    getkey_double("NLOPT_etol", input_vec2d,myinput.NLOPT_etol);
}//end input_nlopt

int input_geom_bond(vector<vector<string> >& input_vec2d){
//    const char* finp="input.inp";
//    vector<vector<string> > input_vec2d;
//    readinp(finp,input_vec2d);
//读入atom块
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
}//input_geom_bond 结束

int input_geom_3d(vector<vector<string> >& input_vec2d){

//    const char* finp="input.inp";
//    vector<vector<string> > input_vec2d;
//    readinp(finp,input_vec2d);
//读入atom块
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
//
vector<int> geom_r;
get_geom_xngro("geom_r",myinput.geom_n,input_vec2d,geom_r);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_r[i]=geom_r[i];
vector<int> geom_theta;
get_geom_xngro("geom_theta",myinput.geom_n,input_vec2d,geom_theta);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_theta[i]=geom_theta[i];
vector<int> geom_phi;
get_geom_xngro("geom_phi",myinput.geom_n,input_vec2d,geom_phi);
for(int i=0;i<myinput.geom_n;i++) myinput.geom_phi[i]=geom_phi[i];

//cout<<myinput.geom_phi[0]<<myinput.geom_phi[1]<<myinput.geom_phi[2]<<endl;

}//input_geom_3d 结束



//不同任务的输入文件
int input_general(vector<vector<string> >& input_vec2d){
//    const char* finp="input.inp";
//    vector<vector<string> > input_vec2d;
//    readinp(finp,input_vec2d);
//cout<<myinput.geom_atom[4][0]<<myinput.geom_atom[4][1]<<myinput.geom_atom[4][2];
//读入pot块map<string,int> potmap;
for(int i=0;i<50;i++)cout<<"#";
cout<<"In input_general"<<endl;
int npot=0;
getkey_int("npot", input_vec2d,npot);
getpot(npot,input_vec2d,myinput.potmap);
//精度参数
getkey_int("ip_num", input_vec2d,myinput.ip_num );//插值点数
myinput.nstr_ALG="COBYLA";//设置默认内层优化算法
getkey_str("nstr_ALG",input_vec2d,myinput.nstr_ALG);//内层优化算法
getkey_int("nstr_num", input_vec2d,myinput.nstr_num);
//调用命令定义
#ifdef MINGW
memcpy(myinput.feffbat,"feffwin.Bat",CHAR_LEN);
memcpy(myinput.feffrun,"feffwin.Bat 1>1.log 2>2.log",CHAR_LEN);
memcpy(myinput.fdmneswin,"fdmnes_win32.exe 1>1.log 2>2.log",CHAR_LEN);
#else
memcpy(myinput.feffbat,"feff.sh",CHAR_LEN);
memcpy(myinput.feffrun,"feff.sh >& 1.log",CHAR_LEN);
memcpy(myinput.fdmneswin,"fdmnes >& 1.log ",CHAR_LEN);
#endif // MINGW
}

//各种不同类型的输入函数
int input_feffprm(vector<vector<string> >& input_vec2d){
//    const char* finp="input.inp";
//    vector<vector<string> > input_vec2d;
//    readinp(finp,input_vec2d);
//输入用临时变量
cout<<"In input_feffprm"<<endl;
    string strin;
    char*  cstrin;
    double ini,low,up;
    double din;
    int    iin;
//
    getkey_str("indir",input_vec2d,strin);
    strcpy(myinput.indir,strin.c_str());
//pdbname
    getkey_str("pdbname",input_vec2d,strin);
    cstrin=new char[CHAR_LEN];//initial char*
    strcpy(cstrin,strin.c_str());
//    const char* ccstr;
//    ccstr=strin.c_str();
//    strcpy(cstrin,ccstr);
    constr(myinput.indir,cstrin,myinput.pdbname);
//    cout<<myinput.pdbname;
 //   std::cout<<myinput.pdbname<<std::endl;
//prmname
    getkey_str("prmname",input_vec2d,strin);
    strcpy(cstrin,strin.c_str());
    constr(myinput.indir,cstrin,myinput.prmname);
    //cout<<myinput.prmname;
//
    getkey_double3("es",input_vec2d,myinput.ini_es,myinput.lo_es,myinput.up_es);
    getkey_double3("ra",input_vec2d,myinput.ini_ra,myinput.lo_ra,myinput.up_ra);
}//end feffprm

int input_ff2x(vector<vector<string> >& input_vec2d){
    //模式5 展宽额外读取
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
}//end input_ff2x

int input_feffdiff(vector<vector<string> >& input_vec2d){


//    const char* finp="input.inp";
//    vector<vector<string> > input_vec2d;
//    readinp(finp,input_vec2d);

    //dis_vec2d(input_vec2d);
    string item;
    string strin;
    char*  cstrin;
    cstrin=new char[CHAR_LEN];//initial char*
    double ini,low,up;
    double din;
    int    iin;
    getkey_str("indir",input_vec2d,strin);
    strcpy(myinput.indir,strin.c_str());
//dis_vec2d(input_vec2d);//显示二维string数组
    getkey_str("diff_exp_char",input_vec2d,strin);
    strcpy(cstrin,strin.c_str());
    constr(myinput.indir,cstrin,myinput.diff_exp_char);

//dis_vec2d(input_vec2d);//显示二维string数组
//cout<<"Here"<<endl;
    getkey_str("pdbname",input_vec2d,strin);
    //getkey_str("pdbname",input_vec2d,strin);
    strcpy(cstrin,strin.c_str());
//    const char* ccstr;
//    ccstr=strin.c_str();
//    strcpy(cstrin,ccstr);
constr(myinput.indir,cstrin,myinput.pdbname);
//  std::cout<<myinput.pdbname<<std::endl;


 // std::cout<<myinput.diff_exp_char<<std::endl;
    getkey_str("prm_the_char",input_vec2d,strin);
//strin="Feprm.feff";
    strcpy(cstrin,strin.c_str());
    constr(myinput.indir,cstrin,myinput.prm_the_char);
// std::cout<<myinput.prm_the_char<<std::endl;

    getkey_double3("ex_diff",input_vec2d,myinput.ini_ex_diff,myinput.lo_ex_diff,myinput.up_ex_diff);
  //  std::cout<<"ex_diff"<<myinput.ini_ex_diff<<" "<<myinput.lo_ex_diff<<"  "<<myinput.up_ex_diff<<std::endl;

    getkey_double3("es_diff",input_vec2d,myinput.ini_es_diff,myinput.lo_es_diff,myinput.up_es_diff);
   // std::cout<<"es_diff"<<myinput.ini_es_diff<<" "<<myinput.lo_es_diff<<"  "<<myinput.up_es_diff<<std::endl;
//注意此处没有读入ra_prm的话会使激发态谱输出为0 补入记录12-7
    getkey_double("es_prm", input_vec2d,myinput.es_prm );
    getkey_double("ra_prm", input_vec2d,myinput.ra_prm );
//    getkey_int("nstr_num", input_vec2d,myinput.nstr_num);

    //diff_exp_char Codpabpydiff.dat

}


int input_fdmprm(vector<vector<string> >&input_vec2d){
    cout<<"In input_fdmprm"<<endl;
    string strin;
    char*  cstrin;
    double ini,low,up;
    double din;
    int    iin;

    getkey_int("nopt_conv", input_vec2d,myinput.nopt_conv);//卷积参数拟合的次数
    myinput.n_conv=5;
    getkey_int("n_conv", input_vec2d,myinput.n_conv);//卷积参数个数一般为5

    get_vector_double(myinput.n_conv,"fdm_up_conv",input_vec2d,myinput.fdm_up_conv);
    get_vector_double(myinput.n_conv,"fdm_low_conv",input_vec2d,myinput.fdm_low_conv);
    get_vector_double(myinput.n_conv,"fdm_x0_conv",input_vec2d,myinput.fdm_x0_conv);
}

//结束不同runtype输入函数




//int print_atom1(double* atom){
//    cout<<atom[0]<<" "<<atom[1]<<" "<<atom[2]<<endl;
//}//
//处理输入文件的函数
int getkey_str(const string a,  vector<vector<string> >& input_vec2d,string& key_str){
    int find=0;//是否找到字符串的变量
    //vector<vector<string> > input_vec2d=input_vec2d_in;
if(myinput.display==1) cout<<"we are searing "<<a<<endl;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==a){
               // if(myinput.display==1)cout<<*j<<" is found"<<endl;
                string ss=*(j+1);
                key_str=ss;
               // if(myinput.display==1)cout<<a<<" is "<<ss<<endl;
                find=1;
                if(myinput.display==1) cout<<a<<" is "<<key_str<<endl;
        }//end if
        }//end for i
       //if(myinput.display==0 && find==0) cout<<a<<" not found"<<endl;
//dis_vec2d(input_vec2d);//里面加一个这个就没问题？
//input_vec2d.clear();
//vector<vector<string> >(input_vec2d).swap(input_vec2d);
//input_vec2d.~vector();

    return 0;
}//end getkey_str


int getkey_double(const string a,  vector<vector<string> >& input_vec2d,double& key_double){
    if(myinput.display==1)cout<<"We are searing "<<a<<endl;
    //vector<vector<string> > input_vec2d=input_vec2d_in;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==a){
                string ss=*(j+1);
                key_double=atof(ss.c_str());
                if(myinput.display==1) cout<<a<<" is "<<key_double<<endl;
        }//end if
        }//end for i
    return 0;
}//end getkey_double
int getkey_int(const string a,  vector<vector<string> >& input_vec2d,int& key_int){
    if(a!="display"&&myinput.display==1) cout<<"We are searching "<<a<<endl;
    //vector<vector<string> > input_vec2d=input_vec2d_in;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==a){
                string ss=*(j+1);

                key_int=atoi(ss.c_str());
                if(a!="display"&&myinput.display==1)cout<<a<<" is "<<key_int<<endl;
        }//end if
        }//end for i
    return 0;
}//end getkey_str
//向量读取三个数值ini lowbond upbond的函数
int getkey_double3(const string a,  vector<vector<string> >& input_vec2d ,double& key_1,double & key_2,double & key_3){
        if(myinput.display==1) cout<<"We are searching "<<a<<endl;
        //vector<vector<string> > input_vec2d=input_vec2d_in;
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
        }//end if
        }//end for i
    return 0;
}//end getkey_str
//读取n长向量的函数
int get_vector_double(int n,const string vecname, vector<vector<string> >& input_vec2d,vector<double>& vec_out){
    //cout<<"npot"<<npot<<endl;
        //vector<vector<string> > input_vec2d=input_vec2d_in;
        if(myinput.display==1)cout<<"We are searching "<<vecname<<endl;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==vecname)iblock=i;
        }//end for i
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
}//end get_vector
//
int getpot(int npot, vector<vector<string> >& input_vec2d,map<string,int>& potmap){
    if(myinput.display==1) cout<<"We are searing pot"<<endl;
    //cout<<"npot"<<npot<<endl;
        //vector<vector<string> > input_vec2d=input_vec2d_in;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j=="pot")iblock=i+1;
        }//end for i
        for(int ip=0;ip<npot;++ip){
            vector<string> item=*(iblock+ip);
            string ele=item[0];
            int pot=atoi(item[1].c_str());
            potmap[ele]=pot;
            if(myinput.display==1)cout<<ele<<" "<<potmap[ele]<<endl;
        }

    return 0;
}//end getpot
//关于几何改变模式
int getatom(int nchange, vector<vector<string> >& input_vec2d,vector<vector<double> >& atom_vec2d){
    //cout<<"nchange"<<nchange<<endl;
       // vector<vector<string> > input_vec2d=input_vec2d_in;
       //dis_vec2d(input_vec2d);
       if(myinput.display==1) cout<<"We are searching representative atom "<<endl;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();

        if(*j=="atom"){iblock=i+1;}
        //if(*j=="es_diff"){iblock=i+1;cout<<"find iblock"<<endl;}
        }//end for i
        //*iblock;
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
}//end getatom

int get_geom_group(int geom_n, vector<vector<string> >& input_vec2d,vector<char>& geom_group){
    if(myinput.display==1) cout<<"We are searching atom group char"<<endl;
    //cout<<"npot"<<npot<<endl;
        //vector<vector<string> > input_vec2d=input_vec2d_in;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j=="geom_group")iblock=i;
        }//end for i
        vector<string> item=*iblock;
        for(int ip=1;ip<geom_n+1;++ip){

            char group_char=item[ip][0];
            geom_group.push_back(group_char);
            if(myinput.display==1) cout<<" "<<group_char;
        }
        if(myinput.display==1) cout<<endl;

    return 0;
}//end get_geom_group

int get_geom_xi(int geom_n, vector<vector<string> >& input_vec2d,vector<int>& geom_xi){
    if(myinput.display==1)cout<<"We are searching xi of group"<<endl;
    //cout<<"npot"<<npot<<endl;
        //vector<vector<string> > input_vec2d=input_vec2d_in;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j=="geom_xi")iblock=i;
        }//end for i
        vector<string> item=*iblock;
        for(int ip=1;ip<geom_n+1;++ip){

            int xi=atoi(item[ip].c_str());
            geom_xi.push_back(xi);
            if(myinput.display==1)cout<<" "<<xi;
        }
        if(myinput.display==1)cout<<endl;

    return 0;
}//end get_geom_grou
//读取r theta phi 3d变换时使用
int get_geom_xngro(const string key_char,int geom_n, vector<vector<string> >& input_vec2d,vector<int>& geom_xi){
    if(myinput.display==1)cout<<"We are searching xi of "<<key_char<<endl;
    //读取 r theta phi的数据
        //vector<vector<string> > input_vec2d=input_vec2d_in;
        vector<vector<string> >::iterator iblock;
        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
        vector<string> vv=*i;
        vector<string>::iterator j=vv.begin();
        if(*j==key_char)iblock=i;
        }//end for i
        vector<string> item=*iblock;
        for(int ip=1;ip<geom_n+1;++ip){

            int xi=atoi(item[ip].c_str());
            geom_xi.push_back(xi);
            if(myinput.display==1)cout<<" "<<xi;
        }
        if(myinput.display==1)cout<<endl;

    return 0;
}//end get_geom_grou











//双层vector内 读取和显示函数
int readinp(const char* finp,vector<vector<string> >& input_vec2d){
    ifstream fin(finp);
    string s;
    while( getline(fin,s) )//读到两层vector里面去
    {   vector<string> oneline;
        if(s.length()==0) continue;//跳过空白行
       // cout << "Read from file: " << s << endl;
        boost::char_separator<char> mysep("  ");//防止标点分割字符串
        boost::tokenizer< boost::char_separator<char> > mytok(s,mysep);//mytok 结果
        int i_tok=0;
        vector<string> mytok_vec;

        //遍历结果mytok
        for(boost::tokenizer<boost::char_separator<char> >::iterator i=mytok.begin();i!=mytok.end();++i)
        {i_tok=i_tok+1;
        mytok_vec.push_back(*i);
        };
        //
        if(mytok_vec[0][0]=='#')continue;//跳过分割后第一个字符串中第一个字符为#的
        input_vec2d.push_back(mytok_vec);
    }//while
return 0;
}

int dis_vec2d(vector<vector<string> > input_vec2d){
    int idx=0;
for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
    vector<string> vv=*i;
    std::cout<<"line "<<idx<<" ";
    for(vector<string>::iterator j=vv.begin();j!=vv.end();++j){
            std::cout<<*j<<" ";

    }//for j
    std::cout<<std::endl;
    ++idx;
}//for i
return 0;
}//dis_vec2d 结束

//废弃读取函数
//int getblock_atom(const char* a, vector<vector<string> > input_vec2d,double& key_1,double & key_2,double & key_3){
//        for(vector<vector<string> >::iterator i=input_vec2d.begin();i!=input_vec2d.end();++i){
//        vector<string> vv=*i;
//        vector<string>::iterator j=vv.begin();
//        if(*j==a){
//                string ss1=*(j+1);
//                key_1=atof(ss1.c_str());
//                string ss2=*(j+2);
//                key_2=atof(ss2.c_str());
//                string ss3=*(j+3);
//                key_3=atof(ss3.c_str());
//        }//end if
//        }//end for i
//    return 0;
//}//end getblock_atom
//void setinput(){//在main头部运行
//
//memcpy(myinput.indir,"inmd/",CHAR_LEN);//设置输入文件夹
////基态拟合参数
////memcpy(myinput.pdbname,"Codpabpy.pdb",CHAR_LEN);
////memcpy(myinput.pdbname,"Coprm6diff5.pdb",CHAR_LEN);
////memcpy(myinput.pdbname,"Coprmfit6.pdb",CHAR_LEN);
////memcpy(myinput.pdbname,"feffprm1.pdb",CHAR_LEN);
////memcpy(myinput.pdbname,"feff1_diffcor5.pdb",CHAR_LEN);
////
////memcpy(myinput.pdbname,"CoIcor6spin3.pdb",CHAR_LEN);
//constr(myinput.indir,"ex-CoImd.pdb",myinput.pdbname);
//
////memcpy(myinput.prmname,"Codpabpyprm.dat",CHAR_LEN);//基态实验谱文件
//constr(myinput.indir,"Codpabpyprm.dat",myinput.prmname);
//
//
//myinput.up_es=30;myinput.lo_es=-30;
//myinput.ini_es=0;
//
//myinput.up_ra=1.5;myinput.lo_ra=0.7;
//myinput.ini_ra=1;
//
//
//
////差分谱拟合相关
////myinput.up_es_diff=5;myinput.lo_es_diff=-5;myinput.ini_es_diff=0;
//myinput.up_es_diff=0;myinput.lo_es_diff=0;myinput.ini_es_diff=0;//固定了能量相对位移
//myinput.up_ex_diff=0.99;myinput.lo_ex_diff=0.01;myinput.ini_ex_diff=0.03;
////memcpy(myinput.diff_exp_char,"Codpabpydiff.dat",CHAR_LEN);//基态实验谱文件
//constr(myinput.indir,"Codpabpydiff.dat",myinput.diff_exp_char);
////memcpy(myinput.prm_the_char,"feff1.prm",CHAR_LEN);//基态理论谱文件
////memcpy(myinput.prm_the_char,"dft.feff",CHAR_LEN);//基态理论谱文件
//constr(myinput.indir,"fms8.feff",myinput.prm_the_char);
//
////-1.37763_34.1535
////myinput.es_prm=-1.37763;myinput.ra_prm=34.1535;
//
////-2.57252_0.931177
////myinput.es_prm=-2.57252;myinput.ra_prm=0.931177;
////dft -3.91154_0.934409
////myinput.es_prm=-3.91154;myinput.ra_prm=0.934409;
////md_es_ratio_-1.62241_0.925563
//myinput.es_prm=-1.62241;myinput.ra_prm=0.925563;
////精度控制参数
//myinput.ip_num=1000;//插值次数
//myinput.nstr_num=1000;// ra es拟合的次数
//
////原子势符号定义
////fdmnes
////myinput.potmap["Co"]=1;
////myinput.potmap["N"]=2;
////myinput.potmap["C"]=3;
////myinput.potmap["H"]=4;
////myinput.potmap["O"]=5;
////feff pot
//myinput.potmap["Co"]=0;
//myinput.potmap["N"]=1;
//myinput.potmap["C"]=2;
//myinput.potmap["H"]=3;
//myinput.potmap["O"]=4;
//
//
////初始化配体形式变动需要的原子坐标
//int n_rep_atom=6;//代表原子数
//myinput.atom_array=(double**)malloc(n_rep_atom*sizeof(double*));
//for(int i=0;i<n_rep_atom;i++){
//    myinput.atom_array[i]=(double*)malloc(3*sizeof(double));
//}
//double a3[3]={-1.131,-1.579,0.287};
//double a4[3]={1.361,-1.349,-0.410};
//myinput.atom_array[0]=a3;myinput.atom_array[1]=a4;
//double a2[3]={0.477,0.053,1.912};
//double a6[3]={-1.357,1.312,0.559};
//myinput.atom_array[2]=a2;myinput.atom_array[3]=a6;
//double a5[3]={1.190,1.513,-0.425};
//double a7[3]={-0.458,0.156,-1.918};
//myinput.atom_array[4]=a5;myinput.atom_array[5]=a7;
//
//
//
////卷积相关参数
//myinput.n_conv=5;
//myinput.up_conv=(double*)malloc(myinput.n_conv*sizeof(double));
//myinput.low_conv=(double*)malloc(myinput.n_conv*sizeof(double));
//myinput.x0_conv=(double*)malloc(myinput.n_conv*sizeof(double));
//
////卷积可调参数
////efermi gamma_hole gamma_max ecent elarg
////myinput.nopt_conv=60;//卷积拟合的次数
////double up_conv[5]={  10,1.33,20,80,50};
////double low_conv[5]={-10,1.33, 0,10, 5};
////double x0_conv[5]={  10,1.33,20,80,50};
//
////固定基态卷积 尤其费米面
//myinput.nopt_conv=2;
////efermi gamma_hole gamma_max ecent elarg
////4.6341,1.33,15.6017,18.4689,49.9004
//double up_conv[5]={4.6341,1.33,15.6017,18.4689,49.9004};
//double low_conv[5]={4.6341,1.33,15.6017,18.4689,49.9004};
//double x0_conv[5]={4.6341,1.33,15.6017,18.4689,49.9004};
//
////
//for(int i=0;i<myinput.n_conv;i++){
//    myinput.up_conv[i]=up_conv[i];
//    myinput.low_conv[i]=low_conv[i];
//    myinput.x0_conv[i]=x0_conv[i];
//}
////结束卷积相关参数
//
////几何结构变化相关
//myinput.geom_model=1;//模式1 刚体 模式2 分层 模式3 刚体加角度
//int geom_n=5;
//myinput.geom_n=5;
////atom 初始化
////double atom2[3]={0.859,-1.195,1.404};
////double atom3[3]={0.275,-1.726,-1.219};
////double atom4[3]={-1.949,-0.467,0.667};
////double atom5[3]={1.918,0.774,-0.430};
////double atom6[3]={-0.053,1.489,1.507};
//
////dft 代表原子
//double atom2[3]={-1.116,0.525,1.527};
//double atom3[3]={-1.831,-0.148,-0.845 };
//double atom4[3]={0.878,-1.871,0.029};
//double atom5[3]={1.154,1.559,-0.685};
//double atom6[3]={1.532,0.225,1.656};
//
//
//
//
//
//
//
//double geom_atom[5][3]={{0.859,-1.195,1.404},{0.275,-1.726,-1.219},{-1.949,-0.467,0.667},{1.918,0.774,-0.430},{-0.053,1.489,1.507}};
//for(int i=0;i<5;i++){
//    for(int j=0;j<3;j++)myinput.geom_atom[i][j]=geom_atom[i][j];
//}
////group初始化
//char group[5]={'B','C','D','E','F'};
//for(int i=0;i<geom_n;i++)myinput.geom_group[i]=group[i];
////优化参数x设置
//int xi[5]={0,1,2,2,2};
//for(int i=0;i<geom_n;i++)myinput.geom_xi[i]=xi[i];
//
//
////调用命令定义
//#ifdef MINGW
//memcpy(myinput.feffbat,"feffwin.Bat",CHAR_LEN);
//memcpy(myinput.feffrun,"feffwin.Bat 1>1.log 2>2.log",CHAR_LEN);
//memcpy(myinput.fdmneswin,"fdmnes_win32.exe 1>1.log 2>2.log",CHAR_LEN);
//#else
//memcpy(myinput.feffbat,"feff.sh",CHAR_LEN);
//memcpy(myinput.feffrun,"feff.sh >& 1.log",CHAR_LEN);
//#endif // MINGW
//
//}//end  setinput

#endif // GLO
