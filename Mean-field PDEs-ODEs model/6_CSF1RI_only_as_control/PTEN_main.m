function PTEN_main(pten)



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
global dmax1 dmax2 K_CSF1RI K_Gal9I K_A B_CSF1RI A_0 B_3
global select_AKT_I dmax3 select_EGFR_I dmax4 select_IGF1R_I dmax5
global select_ERK_I dmax6 select_GSK3_I dmax7 select_IRF1_I dmax8

Parameters_nondimensional
PTEN=pten
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

select_CSF1R_I=0;%%当选择情况自适应情况5的时候，需要用pdepe3来解
select_Gal9_I=0;%2
select_AKT_I=0;%3
select_EGFR_I=0;%4
select_IGF1R_I=0;%5
select_ERK_I=0;%6
select_GSK3_I=0;
select_IRF1_I=0;%8


dmax1=2;%CSF1
dmax2=1;%Gal9
dmax3=1;%AKT
dmax4=3;%EGFR
dmax5=3;%IGF1R
dmax6=3;%ERK
dmax7=3;%GSK
dmax8=3;%IRF1


% global b1 b2 b3 b4
% global t51 t52 t53 t54 t55 t56 t57 t58
% b1=0;b2=1;b3=1;b4=0;
% t51=0*delt_t;t52=30*delt_t;t53=31*delt_t;t54=300*delt_t;t55=100*delt_t;t56=149*delt_t;t57=150*delt_t;t58=201*delt_t;
%
% global b t5 b00 c00 fin1 fin2 tumor_mean_min
% b=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
% b=b*1;
% t50=10000*ones(40,1);
% t5=t50*delt_t;
% t5(1)=0;
% b00=1;
% c00=1;
% fin1=0.3;
% fin2=0.05;
% tumor_mean_min=0;


options=odeset('reltol',1e-20);%精度函数不能少
m = 2;%用matlab解需要设置的参数
sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);% odex4ode：计算希尔方程；pdex4ic：除药物方程外的初始值
u1 = sol(:,:,1);%C_T
u2 = sol(:,:,2);%M1
u3 = sol(:,:,3);%M2
u4 = sol(:,:,4);%CSF1
u5 = sol(:,:,5);%Gal9
u6 = sol(:,:,6);%EGF
u7 = sol(:,:,7);%IGF1
u8 = sol(:,:,8);%EGFR
u9 = sol(:,:,9);%ERK
u10 = sol(:,:,10);%IGF1R
u11 = sol(:,:,11);%AKT
u12 = sol(:,:,12);%GSK
u13 = sol(:,:,13);%IRF1
u14 = sol(:,:,14);%Gal9_细胞内信号通路
C_T_mean=sum(u1(:,:),2)/100;
% M1_mean=sum(u2(:,:),2)/100;
% M2_mean=sum(u3(:,:),2)/100;
% t6=t5/delt_t;
TimeLength=200;
if PTEN==1
    save('.\total_cell\WT\0.mat','C_T_mean');
else
    save('.\total_cell\Loss\0.mat','C_T_mean');%9，16，22，27，31，6-34，36
end

% figure
% pcolor(x,0:1:60,u1(1:61,:))
% title('tumor','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% clim([0 1]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u2(1:61,:))
% title('M1','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.25]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u3(1:61,:))
% title('M2','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.3]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u4(1:61,:))
% title('CSF1','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.5]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u5(1:61,:))
% title('Gal-9','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.5]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u6(1:61,:))
% title('EGF','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.5]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u7(1:61,:))
% title('IGF1','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.25]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u8(1:61,:))
% title('EGFR','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.7]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% 
% figure
% pcolor(x,0:1:60,u9(1:61,:))
% title('ERK','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.7]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u10(1:61,:))
% title('IGF1R','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.14]);
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');
% 
% figure
% pcolor(x,0:1:60,u11(1:61,:))
% title('AKT','Fontweight','Bold','FontSize',18,'FontName','Arial');
% colorbar;
% shading interp;
% % clim([0 0.8])
% set(gca,'ytick',0:10:60); set(gca,'yticklabel',0:10:60,'FontName','Arial');
% set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
% xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
% set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',16,'FontName','Arial');

tmp=0.027;
tmp1=0.04*19;
tmp2=0.300;

% if PTEN==1
%     save('.\参数敏感\WT\1_C.mat','C_T_mean');
% else
%     save('.\参数敏感\Loss\1_C.mat','C_T_mean');
% end

    % figure
    % plot(0:1:200,C_T_mean(1:201),'-','LineWidth',3,'Color',[0.85,0.33,0.10]);hold on;
    % % plot(x_day,y_volume,'s','LineWidth',5)
    % % title('PTEN\_WT IGF1R-I','FontWeight','Bold','FontSize',18)
    % set(gca,'xtick',(0:50:200));
    % ylim([0 1]);
    % % set(gca,'xticklabel',0:5:45,'FontName','Arial');
    % xlabel('Days','FontWeight','Bold','FontSize',18);
    % ylabel('Average Cancer Cell Density','FontWeight','Bold','FontSize',18);
    % set(gca,'FontWeight','Bold','FontSize',16);
    % hold on


% x_day=[8 22];
% y_volume=[0.01 0.01*19];







