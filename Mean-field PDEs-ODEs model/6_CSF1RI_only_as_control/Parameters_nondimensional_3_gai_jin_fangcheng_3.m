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
D_drug=2.0*10^(-6);  %药物的扩散速率
eta_drug=1.0*10^(-17); %药物的降解速率%s-1 %=1 day-1 (Cristini, Cancer Research, 2009.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scaling parameters
t_s=D_IGF1/(d_IGF1*D_T);  %% 92.5926 day
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
B_2=t_s*S_A;
B_3=d_EGF/d_IGF1;
epsilon=D_T/D_IGF1;
D_d=D_drug/D_T;
eta_d=eta_drug*t_s;
r_T=r_T_0*(D_IGF1/d_IGF1/D_T);
r_T=r_T/100;
d_T=d_T_0*(D_IGF1/d_IGF1/D_T);
d_T=d_T/100;
%调试中
b_EGF=1.5;
% a_r=10;%本来应该是10的
a_ERK=20;%本来应该是10的;%
a_AKT=a_ERK*0.4;%本来应该是10的
k_r=2.5;%最终模型用的这个参数
k_r_0=k_r/C_M2_max;%Word的调整参数 1.2500*10^(-8)
k_d=1.3*a_ERK*r_T/d_T/2;%最终模型用的这个参数
% k_d=500;%最终模型用的这个参数
k_d_0=k_d/C_M1_max;%Word的调整参数 7.5000*10^(-7)
k_M=3;%最终模型用的这个参数
k_M_0=k_M/(D_IGF1/d_IGF1/D_T);%Word的调整参数  1.2500e-07\
a_CSF1=k_M*2.5;%最终模型用的这个参数%原来是1.5
a_CSF1_0=a_CSF1/(D_IGF1/d_IGF1/D_T);%Word的调整参数  1.8750e-07
a_A=k_M*2.5;%最终模型用的这个参数
a_A_0=a_A/(D_IGF1/d_IGF1/D_T);%Word的调整参数   2.5000e-07
K1=0.5;%根据最大IGF1设计 最终模型用的这个参数
K1_0=K1/IGF1_s;%Word的调整参数 1.6276e+11
K210=0.3;%根据最大CSF1设计 最终模型用的这个参数     %原来是10
K210_0=K210/CSF1_s;%Word的调整参数 1.5625e+10
K220=10;%根据药物设计 最终模型用的这个参数
K220_0=K220/CSF1_s;%Word的调整参数 3.9062e+10
K3=3;%根据药物的积累设计 最终模型用的这个参数
K5=0.5;
K6=0.5;
%初始化0时刻的耐药细胞因子的量
K4=1;
A_0=0; %最终模型用的这个参数 normalized, IL4 level in Veh tumor cells, ref to Science,2015
x_0=sqrt(D_IGF1/d_IGF1);

load('D:\Matlab代码存放\2018胶质瘤\总程序\canshu114-24hours-ode15s.mat','p');
Ptest1=p(1,1:22);
V1_0=Ptest1(1);
a_0=Ptest1(2);
V21=Ptest1(3);
V22=Ptest1(4);
V3_0=Ptest1(5);
V4_0=Ptest1(6);
d1_0=Ptest1(7);
d2_0=Ptest1(8);
d3_0=Ptest1(9);
d4_0=Ptest1(10);
K11_0=Ptest1(11);
K12=Ptest1(12);
K21=Ptest1(13);
K22=Ptest1(14);
K23=Ptest1(15);
K31_0=Ptest1(16);
K32=Ptest1(17);
K41=Ptest1(18);
K42=Ptest1(19);
K13=Ptest1(20);
K33=Ptest1(21);
n=floor(Ptest1(22));
V1=V1_0*t_s;
a=a_0*t_s;
V3=V3_0*t_s;
V4=V4_0*t_s;
% K11=K11_0/EGF_s;
K11=K11_0;
% K31=K31_0/IGF1_s;
K31=K31_0;
d1=d1_0*t_s;
d2=d2_0*t_s;
d3=d3_0*t_s;
d4=d4_0*t_s;