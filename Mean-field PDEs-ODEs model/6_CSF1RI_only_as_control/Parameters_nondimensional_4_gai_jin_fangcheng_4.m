%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Original parameters
D_T=1.0*10^(-9); % cm2/s;  %6.0 um2/min   %肿瘤细胞的扩散率
D_M1=1.0*10^(-9); % cm2/s;  %6.0 um2/min   %M1巨噬细胞的扩散率
D_M2=1.0*10^(-9); % cm2/s;  %6.0 um2/min   %M2巨噬细胞的扩散率
r_T_0=2.41*10^(-5);  %肿瘤细胞的基础增值率
d_T_0=2.8660*1e-7;  %肿瘤细胞的基础死亡率
a_T_EGF=1.7*10^(-3); %表示肿瘤细胞对于EGF的趋化系数
a_T_IGF1=1.7*10^(-3); %表示肿瘤细胞对于IGF1的趋化系数
a_M1=a_T_IGF1;%M1巨噬细胞对于CSF1的趋化系数
a_M2=a_T_IGF1;%M2巨噬细胞对于CSF1的趋化系数
D_CSF1=2.0*10^(-6);  % cm2s-1            %CSF1的扩散率
D_EGF=2.0*10^(-6);                      %EGF的扩散率
D_IGF1=2.0*10^(-6);                      %IGF1的扩散率
S_CSF1=3.20*10^(-22);  % gs-1cell-1      %肿瘤细胞对于CSF1的分泌率
S_EGF=1.92*10^(-22);  % gs-1cell-1      %M2巨噬细胞对于EGF的分泌率
S_IGF1=1.92*10^(-22);  % gs-1cell-1      %M2巨噬细胞对于IGF1的分泌率
%CSF1的降解率
d_CSF1=5.0*10^(-4); %4.80*10^(-4) ; %%3*10^(-2)/60 ; % s-1 (JTB)   %4.80*10^(-5) %s-1  (BMB paper)
%EGF的降解率
d_EGF=2.5*10^(-4); 
%IGF1的降解率
d_IGF1=2.5*10^(-4); % s-1 (JTB)  % 2.00*10^(-5); %s-1  (BMB)
%CSF1_R_I药物累计产生的耐药细胞因子的速率
S_A=1.2905*10^(-7);
D_drug_old=2.0*10^(-6);  %药物的扩散速率%旧的，后面D_drug是新的
eta_drug=1.0*10^(-17); %药物的降解速率%s-1 %=1 day-1 (Cristini, Cancer Research, 2009.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scaling parameters
t_s=D_IGF1/(d_IGF1*D_T);  %% 92.5926 day
x_0=sqrt(D_IGF1/d_IGF1);  %%无量纲化用到的x
D_drug=D_drug_old*x_0^2;%药物的扩散速率%新的
C_T_max=1.0*10^9;%肿瘤密度的最大值
C_M1_max=2.0*10^8;%M1巨噬细胞密度的最大值
C_M2_max=2.0*10^8;%M2巨噬细胞密度的最大值
CSF1_s=S_CSF1/d_CSF1*C_T_max;
EGF_s=S_EGF/d_EGF*C_M1_max;
IGF1_s=S_IGF1/d_IGF1*C_M1_max;
CSF1_max=1;
CSF1_max_0=CSF1_max*CSF1_s;
EGF_max=1;
EGF_max_0=EGF_max*EGF_s;
IGF1_max=1;
IGF1_max_0=IGF1_max*IGF1_s;

A_T_EGF=a_T_EGF/D_T*S_EGF*C_M1_max/d_EGF;
A_T_IGF1=a_T_IGF1/D_T*S_IGF1*C_M1_max/d_IGF1;
A_M1=a_M1/D_M1*S_CSF1/d_CSF1*C_T_max;
A_M2=a_M2/D_M2*S_CSF1/d_CSF1*C_T_max;
B_1=d_CSF1/d_IGF1;
B_2=d_EGF/d_IGF1;
B_3=t_s*S_A;
epsilon=D_T/D_IGF1;
% D_d=D_drug/D_T;
D_d=D_drug*t_s/x_0^2;
eta_d=eta_drug*t_s;
r_T=r_T_0*(D_IGF1/d_IGF1/D_T);
r_T=r_T/100;
d_T=d_T_0*(D_IGF1/d_IGF1/D_T);
d_T=d_T/100;
%调试中
b_EGF=1.5;
% a_r=10;%本来应该是10的
a_ERK=20;%本来应该是10的;%
a_AKT=a_ERK*0.3;%本来应该是10的
k_r=2.5;%最终模型用的这个参数
k_r_0=k_r/C_M2_max;%Word的调整参数 1.2500*10^(-8)
d_TM1=1.3*a_ERK*r_T/d_T/2;%最终模型用的这个参数
% d_TM1=500;%最终模型用的这个参数
d_TM1_0=d_TM1/C_M1_max;%Word的调整参数 7.5000*10^(-7)
a_M2M1=10;%最终模型用的这个参数
a_M2M1_0=a_M2M1/(D_IGF1/d_IGF1/D_T);%Word的调整参数  1.2500e-07\
a_CSF1=a_M2M1*2.5;%最终模型用的这个参数%原来是1.5
a_CSF1_0=a_CSF1/(D_IGF1/d_IGF1/D_T);%Word的调整参数  1.8750e-07
a_A=a_M2M1*8;%最终模型用的这个参数
a_A_0=a_A/(D_IGF1/d_IGF1/D_T);%Word的调整参数   2.5000e-07
K1=0.5;
K2=0.5;
K31=0.3;%根据最大CSF1设计 最终模型用的这个参数     %原来是10
K31_0=K31/CSF1_s;%Word的调整参数 1.5625e+10
K32=10;%根据药物设计 最终模型用的这个参数
K32_0=K32/CSF1_s;%Word的调整参数 3.9062e+10
K4=30;%根据药物的积累设计 最终模型用的这个参数

A_0=0; %最终模型用的这个参数 normalized, IL4 level in Veh tumor cells, ref to Science,2015

% load('D:\Matlab代码存放\2018胶质瘤\总程序\canshu114-24hours-ode15s.mat','p');
load('canshu1-24hours-ode15s.mat','p');
Ptest1=p(1,1:22);
V7_0=Ptest1(1);
a_0=Ptest1(2);
V81=Ptest1(3);
V82=Ptest1(4);
V9_0=Ptest1(5);
V10_0=Ptest1(6);
d7_0=Ptest1(7);
d8_0=Ptest1(8);
d9_0=Ptest1(9);
d10_0=Ptest1(10);
K71_0=Ptest1(11);
K72=Ptest1(12);
K81=Ptest1(13);
K82=Ptest1(14);
K83=Ptest1(15);
K91_0=Ptest1(16);
K92=Ptest1(17);
K101=Ptest1(18);
K102=Ptest1(19);
K73=Ptest1(20);
K93=Ptest1(21);
n=floor(Ptest1(22));
V7=V7_0*t_s;
a=a_0*t_s;
V9=V9_0*t_s;
V10=V10_0*t_s;
% K71=K71_0/EGF_s;
K71=K71_0;
% K91=K91_0/IGF1_s;
K91=K91_0;
d7=d7_0*t_s;
d8=d8_0*t_s;
d9=d9_0*t_s;
d10=d10_0*t_s;
a_M2M1=a_M2M1*1.4400;
a_CSF1=a_CSF1*5;
a_CSF1_0=a_CSF1/(D_IGF1/d_IGF1/D_T);
% a_ERK=a_ERK*0.5;
a_ERK=a_ERK*0.65*1.1;
% a_AKT=a_AKT/0.3*2;
a_AKT=a_AKT/0.3*2.2;
d_TM1=d_TM1*1.3;
d_TM1_0=d_TM1/C_M1_max;


a_A=2.1*a_A;
a_A_0=a_A/(D_IGF1/d_IGF1/D_T);%Word的调整参数   2.5000e-07