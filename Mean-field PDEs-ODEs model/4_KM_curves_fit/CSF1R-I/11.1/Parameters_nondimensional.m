%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Original parameters
D_T=1.0*10^(-9); % cm2/s;  %6.0 um2/min   %肿瘤细胞的扩散率
D_M1=1.0*10^(-9); % cm2/s;  %6.0 um2/min   %M1巨噬细胞的扩散率
D_M2=1.0*10^(-9); % cm2/s;  %6.0 um2/min   %M2巨噬细胞的扩散率
D_CSF1=2.0*10^(-6);  % cm2s-1            %CSF1的扩散率
D_Gal9=2.0*10^(-6);                      %Gal9的扩散率
D_EGF=2.0*10^(-6);                      %EGF的扩散率
D_IGF1=2.0*10^(-6);%IGF1的扩散率
D_drug=2.0*10^(-6);

r_T_0=2.41*10^(-5)/2/2;  %肿瘤细胞的基础增值率
d_T_0=2.8660*1e-7;  %肿瘤细胞的基础死亡率


alpha_EGF=1.7*10^(-3); %表示肿瘤细胞对于EGF的趋化系数
alpha_IGF1=1.7*10^(-3); %表示肿瘤细胞对于IGF1的趋化系数
alpha_CSF1=1.7*10^(-3);%M1巨噬细胞对于CSF1的趋化系数
alpha_Gal9=1.7*10^(-3);%M1巨噬细胞对于Gal9的趋化系数

S_CSF1=3.20*10^(-22);  % gs-1cell-1      %肿瘤细胞对于CSF1的分泌率
S_Gal9_0=3.20*10^(-22);%Gal9
S_EGF=1.92*10^(-22);  % gs-1cell-1      %M2巨噬细胞对于EGF的分泌率
S_IGF1=1.92*10^(-22);  % gs-1cell-1      %M2巨噬细胞对于IGF1的分泌率

d_CSF1=5.0*10^(-4); %4.80*10^(-4) ; %%3*10^(-2)/60 ; % s-1 (JTB)   %4.80*10^(-5) %s-1  (BMB paper)
d_Gal9=5.0*10^(-4);
d_EGF=2.5*10^(-4);
d_IGF1=2.5*10^(-4); % s-1 (JTB)  % 2.00*10^(-5); %s-1  (BMB)
eta_drug=8.02*10^(-6); %药物的降解速率%s-1 %=1 day-1 (Cristini, Cancer Research, 2009.)
C_T_max=1.0*10^9;%肿瘤密度的最大值
C_M_max=2.0*10^8;%M1,M2巨噬细胞密度的最大值

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以下为估计

CSF1_max=6.4*10^(-10);
EGF_max=1.5360*10^(-10);
IGF1_max=1.5360*10^(-10);

alpha_ERK=14.3;%本来应该是10的;%
alpha_AKT=44;%本来应该是10的
% alpha_AKT=10000;%本来应该是10的

d_TM1_0=1.6397*10^(-5);

beta_CSF1=4.6875*10^(-6);
beta_Gal9=1.75*10^(-5); %Gal9
% beta_M1M2=1.25*10^(-6);
beta_M1M2=0.27*10^(-6);

V_Gal9=2.5*10^(-5);
S_Gal9=S_Gal9_0*V_Gal9;

K1=1.92*10^(-10);%CSF1
K2=0.3*10^(-13);%Gal9_en
K3=0.0001;%Gal9
K4=0.5;%ERK
K5=0.25;%AKT

S_A=1.2905*10^(-7)*3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scaling parameters
t_s=D_IGF1/(d_IGF1*D_T);  %% 92.5926 day
x_0=sqrt(D_IGF1/d_IGF1);  %%无量纲化用到的x

EGF_s=S_EGF/d_EGF*C_M_max;
IGF1_s=S_IGF1/d_IGF1*C_M_max;
CSF1_s=S_CSF1/d_CSF1*C_T_max;
Gal9_s=S_Gal9/d_Gal9*C_T_max;

A_EGF=alpha_EGF/D_T*S_EGF*C_M_max/d_EGF;
A_IGF1=alpha_IGF1/D_T*S_IGF1*C_M_max/d_IGF1;
A_CSF1=alpha_CSF1/D_M1*S_CSF1*C_T_max/d_CSF1;
A_Gal9=alpha_Gal9/D_M1*S_Gal9*C_T_max/d_Gal9;

B_CSF1=D_IGF1/d_IGF1/D_T*beta_CSF1;
B_Gal9=D_IGF1/d_IGF1/D_T*beta_Gal9;
B_M1M2=D_IGF1/d_IGF1/D_T*beta_M1M2;
epsilon=D_T/D_IGF1;
M1=d_CSF1/d_IGF1;
M2=d_Gal9/d_IGF1;
M3=d_EGF/d_IGF1;

D_d=D_drug/D_T;
eta_d=eta_drug*t_s;

CSF1_max=1;
EGF_max=1;
IGF1_max=1*5;

r_T=r_T_0*t_s;
r_T=r_T/100;
d_T=d_T_0*t_s;
d_T=d_T/100;
d_TM1=C_M_max*d_TM1_0;%最终模型用的这个参数

K1=K1/CSF1_s;
K2=K2/Gal9_s/3;

b_EGF=1.5;

k_r=2.5;%最终模型用的这个参数
k_r_0=k_r/C_M_max;%Word的调整参数 1.2500*10^(-8)
a_A=a_M2M1*8;%最终模型用的这个参数
a_A_0=a_A/(D_IGF1/d_IGF1/D_T);%Word的调整参数   2.5000e-07


a_M2M1=10;%最终模型用的这个参数
a_M2M1_0=a_M2M1/(D_IGF1/d_IGF1/D_T);%Word的调整参数  1.2500e-06\
a_CSF1=a_M2M1*2.5;%最终模型用的这个参数%原来是1.5
a_CSF1_0=a_CSF1/(D_IGF1/d_IGF1/D_T);%Word的调整参数  1.8750e-07

A_0=0; %最终模型用的这个参数 normalized, IL4 level in Veh tumor cells, ref to Science,2015

% load('D:\Matlab代码存放\2018胶质瘤\总程序\canshu114-24hours-ode15s.mat','p');
load('pa.mat','Ptest1');
Ptest=Ptest1;
p_ind=[1 2 5 6 7 9 10 11 12 13 14 15 16 17];
for i = p_ind
    Ptest(i)=Ptest(i)*t_s;
end
    
Ptest(18)=Ptest(18);
Ptest(23)=Ptest(23);



