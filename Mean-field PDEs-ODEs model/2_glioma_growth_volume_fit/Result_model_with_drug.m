% addpath('/home/lhf/matlab code/funfun')%add operation path
global eta1    %药物的降解率
global D_d1    %药物的扩散率
global delt_t  %模拟时候的时间间隔
global L      %模拟时候横轴x的大小
global B_3  %这个是需要估计的值,是IL4分泌速率相关的
global A_0   %这个是需要估计的值,是耐药细胞因子初始值
global period1       %真实的的时间周期，由period×delt_t得到
global inter1        %真实的的时间间隔，由inter×delt_t得到
global TimeLength    %总的时间
global Ptest1
global EGF IGF1
global EGFR_I IGF1R_I
global C_T_0 C_M1_0 C_M2_0 CSF1_0 EGF_0 IGF1_0
global ERK AKT
global x_0
global V7 a V81 V82 V9 V10 d7 d8 d9 d10 K71 K72 K81 K82 K83 K91 K92 K101 K102 K73 K93 n
global A_T_EGF A_T_IGF1 A_M1 A_M2 epsilon r_T a_ERK K1 K2 K31 K32 K4 a_AKT  d_T d_TM1 a_M2M1 a_CSF1  a_A K3 B_1 B_2 b_EGF
global CSF1_max EGF_max IGF1_max
global select_CSF1R_I select_EGFR_I select_IGF1R_I 
Parameters_nondimensional_4_gai_jin_fangcheng_4
% load('D:\Matlab代码存放\2018胶质瘤\总程序\canshu114-24hours-ode15s.mat','p');
load('canshu1-24hours-ode15s.mat','p');
Ptest1=p(1,1:22);
TimeLength=200;
L=100*x_0;         %模拟时候坐标轴的总长度
% x_gap=1;                %模拟时候坐标轴的间隔
x = 0:x_0:L;            %肿瘤生长的空间，一维的
% x = 0:1:L;            %肿瘤生长的空间，一维的
delt_t=1/91.5926/24;      %模拟时候的时间间隔
% delt_t=1/91.5926;      %模拟时候的时间间隔
%delt_t是一天，因为3600*24/t_s= 0.0108=1/92.5926
t=0:delt_t:delt_t*24*TimeLength;
% t=0:delt_t:delt_t*100;
% period=20;           %小周期，例如20天用药量6，20天用药量0
period1=period*delt_t*24;  %真正周期
inter=1;             %间隔
inter1=inter*delt_t*24;   %真正间隔
eta1=eta_d;
D_d1=D_d;%药物的扩散率


select_CSF1R_I=1;
select_EGFR_I=3;
select_IGF1R_I=0;



options=odeset('reltol',1e-3);%精度函数不能少
% options=odeset('reltol',1);%精度函数不能少
m = 2;%用matlab解需要设置的参数
sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);