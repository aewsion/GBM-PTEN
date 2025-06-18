function sum_1



addpath('.\funfun');


global eta1    %药物的降解率
global D_d1    %药物的扩散率
global delt_t  %模拟时候的时间间隔
global L  L2    %模拟时候横轴x的大小
global B_3  %这个是需要估计的值,是IL4分泌速率相关的
global A_0   %这个是需要估计的值,是耐药细胞因子初始值
global d_max         %药物边界值，等于3
global period        %这个是用于边界用药时候有一个周期，例如7天一个周期
global inter         %阶梯函数的换化间隔，例如用一天把药物量为0渐变成药物量为6，一天就是这个inter
global period1 period_2      %真实的的时间周期，由period×delt_t得到 period_2是震荡函数的周期
global inter1        %真实的的时间间隔，由inter×delt_t得到
global TimeLength    %总的时间
global m0        %边界条件三设计的药物的量
global Ptest
global EGF IGF1
global EGFR_I IGF1R_I
global C_T_0 C_M1_0 C_M2_0 CSF1_0 EGF_0 IGF1_0 Gal9_0
global ERK AKT
global x_0 %模拟时候坐标轴的总长度 x_0=sqrt(D_IGF1/d_IGF1);
% global V7 a V81 V82 V9 V10 d7 d8 d9 d10 K71 K72 K81 K82 K83 K91 K92 K101 K102 K73 K93 n
global A1 A2 A3 A4 A5 A6 B1 B2 B3 epsilon M1 M2 M3 r_T a_ERK a_AKT d_T d_TM1 K1 K2 K3 K4 K5
global CSF1_max EGF_max IGF1_max dmax1 dmax2 dmax3
global select_CSF1R_I select_EGFR_I select_IGF1R_I a_A
global PTEN b_EGF
Parameters_nondimensional
load('pa.mat','Ptest1')
Ptest=Ptest1;
PTEN=0;
% Conditions_ODES_2;\



TimeLength=200;
L2=100;
L=1/x_0 ;  %模拟时候坐标轴的总长度 x_0=sqrt(D_IGF1/d_IGF1);
x = 0:L/L2:L;            %肿瘤生长的空间，一维的
delt_t=1/91.5926;      %模拟时候的时间间隔，delt_t是一天，因为3600*24/t_s= 0.0108=1/92.5926 t_s=D_IGF1/(d_IGF1*D_T);  %% 92.5926 day
t=0:delt_t:delt_t*TimeLength;
period_2=delt_t*14;
period=7;           %小周期，例如20天用药量6，20天用药量0
period1=period*delt_t;  %真正周期
T=period1*2; %真正的周期；
inter=1;             %间隔 插值d
inter1=inter*delt_t;   %真正间隔
eta1=eta_d;
D_d1=D_d;%药物的扩散率


select_CSF1R_I=1;%%当选择情况自适应情况5的时候，需要用pdepe3来解
select_EGFR_I=0;
select_IGF1R_I=0;

dmax1=1;
dmax2=1;
dmax3=1;

a_A=2.5*a_A;

global b1 b2 b3 b4
global t51 t52 t53 t54 t55 t56 t57 t58
b1=0;b2=1;b3=1;b4=0;
t51=0*delt_t;t52=30*delt_t;t53=31*delt_t;t54=300*delt_t;t55=100*delt_t;t56=149*delt_t;t57=150*delt_t;t58=201*delt_t;

global b t5 b00 c00 fin1 fin2 tumor_mean_min
b=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
b=b*1;
t50=10000*ones(40,1);
t5=t50*delt_t;
t5(1)=0;
b00=1;
c00=1;
fin1=0.3;
fin2=0.05;
tumor_mean_min=0;
options=odeset('reltol',1e-20);%精度函数不能少
m = 2;%用matlab解需要设置的参数
sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);% odex4ode：计算希尔方程；pdex4ic：除药物方程外的初始值
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);
u5 = sol(:,:,5);
u6 = sol(:,:,6);
u7 = sol(:,:,7);
u8 = sol(:,:,8);
u9 = sol(:,:,9);
u10 = sol(:,:,10);
u11 = sol(:,:,11);
u12 = sol(:,:,12);
u13 = sol(:,:,13);
u14 = sol(:,:,14);
C_T_mean=sum(u1(:,:),2)/100;
M1_mean=sum(u2(:,:),2)/100;
M2_mean=sum(u3(:,:),2)/100;
t6=t5/delt_t;
TimeLength=200;


x00=0:1:L2;
t00=0:1:TimeLength;



[a b]=find(u1<0)
figure
pcolor(x,t,u1)
title('tumor','Bold','FontSize',18,'FontName','Arial');
colorbar;
shading interp;
set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
set(gca,'FontSize',16,'FontName','Arial');

figure
pcolor(x,t,u2)
title('M1','FontWeight','Bold','FontSize',18,'FontName','Arial');
colorbar;
shading interp;
set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
set(gca,'FontSize',16,'FontName','Arial');



figure
pcolor(x,t,u3)
title('M2','FontWeight','Bold','FontSize',18,'FontName','Arial');
colorbar;
shading interp;
set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
set(gca,'FontSize',16,'FontName','Arial');



figure

pcolor(x,t,u14)
title('Gal9','Bold','FontSize',18,'FontName','Arial');
colorbar;
shading interp;
set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
set(gca,'FontSize',16,'FontName','Arial');





