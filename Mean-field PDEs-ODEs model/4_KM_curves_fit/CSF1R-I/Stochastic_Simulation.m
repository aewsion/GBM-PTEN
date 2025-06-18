%%%%%% Stochastic Simulatin
%%%%%% Xiaoqiang Sun @ SYSU
%%%%%% 2016,10,15
clear all
close all

tic


addpath('.\funfun');


global x_0 delt_t %\hat(x), \hat(t), x_0=sqrt(D_IGF1/d_IGF1);
global L L2 %模拟时候横轴x的大小
global TimeLength
global Ptest
global EGF IGF1
global A_EGF A_IGF1 A_CSF1 A_Gal9
global B_CSF1 B_Gal9 B_M1M2
global epsilon M1 M2 M3
global CSF1_max EGF_max IGF1_max
global r_T d_T d_TM1
global alpha_ERK alpha_AKT
global K1 K2 K3 K4 K5
global PTEN

global B_IL4  %这个是需要估计的值,是IL4分泌速率相关的
global A_0   %这个是需要估计的值,是耐药细胞因子初始值
global d_max         %药物边界值，等于3
global period        %这个是用于边界用药时候有一个周期，例如7天一个周期
global inter         %阶梯函数的换化间隔，例如用一天把药物量为0渐变成药物量为6，一天就是这个inter
global period1 period_2      %真实的的时间周期，由period×delt_t得到 period_2是震荡函数的周期
global inter1        %真实的的时间间隔，由inter×delt_t得到
global m0        %边界条件三设计的药物的量
global EGFR_I IGF1R_I AKT_I GSK_I

global select_CSF1R_I select_Gal9_I
global b_EGF

global eta_d D_d %药物
global dmax1 dmax2 K_CSF1RI K_Gal9I K_A B_CSF1RI A_0 B_3 S_IGF1 D_T S_A

Parameters_nondimensional
PTEN=1;
AKT_I=0;
GSK_I=0;
A_0=0;

TimeLength=200;
L2=100;
L=2/x_0 ;  %模拟时候坐标轴的总长度 x_0=sqrt(D_IGF1/d_IGF1);
x = 0:L/L2:L;            %肿瘤生长的空间，一维的
delt_t=1/91.5926;      %模拟时候的时间间隔，delt_t是一天，因为3600*24/t_s= 0.0108=1/92.5926 t_s=D_IGF1/(d_IGF1*D_T);  %% 92.5926 day
t=0:delt_t:delt_t*TimeLength;
period_2=delt_t*14;
period=7;           %小周期，例如20天用药量6，20天用药量0
period1=period*delt_t;  %真正周期
T=period1*2; %真正的周期；
inter=1;             %间隔 插值d
inter1=inter*delt_t;   %真正间隔

select_CSF1R_I=1;%PTEN=1
select_Gal9_I=0;%PTEN=0

dmax1=1;
dmax2=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters
Parameter=[B_M1M2 d_TM1 r_T S_A];%B_3==S_A

p = 100;   % Number of points
N = 1;   % Number of dimensions

% mu=[1.05*a_ERK 0.05*a_A];
% sigma=[0.005*a_ERK 0;0 0.03*a_A];
% for i=1:1:1000%取大于0的数
%     bm41=mvnrnd(mu,sigma,p);
%     negative_count=sum(sum(bm41<=0));
%     if negative_count==0
%         break;
%     end
% end
%
% mu=[1.05*a_ERK 2.47*a_A];
% sigma=[0.005*a_ERK 0;0 32*a_A];
% for i=1:1:1000
%     bm42=mvnrnd(mu,sigma,p);
%     negative_count=sum(sum(bm42<=0));
%     if negative_count==0
%         break;
%     end
% end
%
% bm4=zeros(p,2);
% rrr=0.6;
% for i=1:1:p
%     rrr1=rand;
%     bm4(i,:)=bm41(i,:)*heaviside(rrr1-rrr)+bm42(i,:)*heaviside(rrr-rrr1);
% end

mu=[0.8*B_M1M2 0.1*B_CSF1RI];%1
sigma=[0.46*B_M1M2 0;0 0.5*B_CSF1RI];%0.01
for i=1:1:1000%取大于0的数
    bm41=mvnrnd(mu,sigma,p);
    negative_count=sum(sum(bm41<=10^(-5)));
    if negative_count==0
        break;
    end
end

mu=[0.8*B_M1M2 2*B_CSF1RI];%num=1全是这行
sigma=[0.46*B_M1M2 0;0 3.5*B_CSF1RI];

for i=1:1:1000%取大于0的数
    bm42=mvnrnd(mu,sigma,p);
    negative_count=sum(sum(bm41<=10^(-5)));
    if negative_count==0
        break;
    end
end

bm4=bm41;
num_1=0.6;
for i=1:1:p
    rand_1=rand;
    bm4(i,:)=bm41(i,:)*heaviside(rand_1-num_1)+bm42(i,:)*heaviside(num_1-rand_1);%200个病人每个都B_M1M2服从正态，100个的S_A是敏感峰，100个S_A是耐药峰
end
% % % % % bm4=(1-rrr)*bm41+rrr*bm42;
r_RG=zeros(1,p);
r_K=zeros(1,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters
C=zeros(p,TimeLength+1);
for i=1:p
    i
    B_M1M2=bm4(i,1);
    B_CSF1RI=bm4(i,2);
    options=odeset('reltol',1e-3);%精度函数不能少
    m = 2;%用matlab解需要设置的参数
    sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);
    u1 = sol(:,:,1);
    C_T_mean=sum(u1(:,:),2)/100;
    C(i,:)=C_T_mean;

    if C(i,TimeLength+1)~=min(C(i,:))
        r_RG(i)=(C(i,TimeLength+1)-min(C(i,:)))/(TimeLength+1-find(C(i,:)==min(C(i,:)),1));
    else
        r_RG(i)=0;
    end
    if C(i,1)~=min(C(i,:))
        r_K(i)=(C(i,1)-min(C(i,:)))/(find(C(i,:)==min(C(i,:)),1));
    else
        r_K(i)=0;
    end
end
if PTEN==0
    save('stochastic_KO.mat','C','r_RG','r_K','bm4');
else
    save('stochastic_WT.mat','C','r_RG','r_K','bm4');
end
survival_fig;






