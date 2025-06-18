function sum_1


% addpath('D:\software\matlabR2022b\toolbox\matlab\funfun'); 
addpath('D:\Matlab代码存放\2018胶质瘤\总程序\funfun'); 
% addpath('D:\Matlab代码存放\2018胶质瘤\总程序\funfun\private'); 

global eta1    %药物的降解率
global D_d1    %药物的扩散率
global delt_t  %模拟时候的时间间隔
global L      %模拟时候横轴x的大小
global B_3  %这个是需要估计的值,是IL4分泌速率相关的
global A_0   %这个是需要估计的值,是耐药细胞因子初始值
global d_max         %药物边界值，等于3
global period        %这个是用于边界用药时候有一个周期，例如7天一个周期
global inter         %阶梯函数的换化间隔，例如用一天把药物量为0渐变成药物量为6，一天就是这个inter
global period1       %真实的的时间周期，由period×delt_t得到
global inter1        %真实的的时间间隔，由inter×delt_t得到
global TimeLength    %总的时间
global m0        %边界条件三设计的药物的量
global Ptest1
global EGF IGF1
global EGFR_I IGF1R_I
global C_T_0 C_M1_0 C_M2_0 CSF1_0 EGF_0 IGF1_0
global ERK AKT
global x_0
global V7 a V81 V82 V9 V10 d7 d8 d9 d10 K71 K72 K81 K82 K83 K91 K92 K101 K102 K73 K93 n
global A_T_EGF A_T_IGF1 A_M1 A_M2 epsilon r_T a_ERK K1 K2 K31 K32 K4 a_AKT  d_T d_TM1 a_M2M1 a_CSF1  a_A K3 B_1 B_2 b_EGF
global CSF1_max EGF_max IGF1_max dmax1 dmax2 dmax3 
global select_CSF1R_I select_EGFR_I select_IGF1R_I 
% Parameters_nondimensional_4_gai_jin_fangcheng_4
Parameters_nondimensional_5_gai_jin_fangcheng_5
% load('D:\Matlab代码存放\2018胶质瘤\总程序\canshu114-24hours-ode15s.mat','p');
load('D:\Matlab代码存放\2018胶质瘤\总程序\canshu1-24hours-ode15s.mat','p');
Ptest1=p(1,1:22);

% Conditions_ODES_2;\



TimeLength=300;
L2=100;
% L=L2*x_0;         %模拟时候坐标轴的总长度
% L=L2 ;  %模拟时候坐标轴的总长度
L=1/x_0 ;  %模拟时候坐标轴的总长度
% x_gap=1;                %模拟时候坐标轴的间隔
% x = 0:x_0:L;            %肿瘤生长的空间，一维的
x = 0:L/100:L;            %肿瘤生长的空间，一维的
delt_t=1/91.5926;      %模拟时候的时间间隔
% delt_t=1/91.5926;      %模拟时候的时间间隔
%delt_t是一天，因为3600*24/t_s= 0.0108=1/92.5926
t=0:delt_t:delt_t*TimeLength;
% t=0:delt_t:delt_t*100;
period=10;           %小周期，例如20天用药量6，20天用药量0%间断用药的周期
period1=period*delt_t;  %真正周期
inter=1;             %间隔
inter1=inter*delt_t;   %真正间隔
eta1=eta_d;
D_d1=D_d;%药物的扩散率


select_CSF1R_I=1;%CSF1RI的开关，1为连续用药，3为间断用药
select_EGFR_I=0;
select_IGF1R_I=0;

dmax1=1;%CSF1RI的用量
% dmax2=0.01;%EGFRI的用量
% dmax3=0.01;%IGF1RI的用量

% global t51 t52 t53
% t51_0=0;
% t52_0=80;
% t53_0=t52_0+inter;
% t51=t51_0*delt_t;
% t52=t52_0*delt_t;
% t53=t53_0*delt_t;



options=odeset('reltol',1e-20);%精度函数不能少
% options=odeset('reltol',1);%精度函数不能少
m = 2;%用matlab解需要设置的参数
sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);
% sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t);
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
C_T_mean=sum(u1(:,:),2)/100; 
% M1_mean=sum(u2(:,:),2)/100; 
% M2_mean=sum(u3(:,:),2)/100; 
% [fff1 fff2]=min(C_T_mean)
% fff3=find(C_T_mean>=0.6);
% fff4=fff3(1,1);
% % save('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\用药CSF1RI一段时间停\CSF1RI_0.1_80_IGF1RI_0.01.mat','u1','C_T_mean');
% load('D:\Matlab代码存放\2018胶质瘤\总程序\新用药组合\yongyaoduibi1_CSF1RI_1.mat','C_T_mean','fff4');

% C_T_meanall(:,9)=C_T_mean;
% fff4_all(1,9)=fff4;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\用药比较\yongyaoduibi1_CSF1RI_0-0.8.mat','C_T_meanall','fff4_all');

% save('D:\Matlab代码存放\2018胶质瘤\总程序\用药比较\yongyaoduibi1_CSF1RI_0.5_xuanze1_EGFRI_0.5_xuanze1.mat','C_T_mean','fff4'); 
% save('D:\Matlab代码存放\2018胶质瘤\总程序\xiangtu311.mat','C_T_mean','M1_mean','M2_mean');  
% save('D:\Matlab代码存放\2018胶质瘤\总程序\xiangtu313.mat','C_T_mean','M1_mean','M2_mean'); 
% B_3=1*B_3;
% threshold=0.55;
% person_0=2;
% C=zeros(person_0+1,24*TimeLength+1);
% C(1,:)=C_T_mean;
% for person=1:1:person_0 
%       B_3=B_3*(1+rand(1));
%       sol = pdepe2(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,options);
%       u1 = sol(:,:,1);
%       C_T_mean=sum(u1(:,:),2)/100; 
%       C(person+1,:)=C_T_mean;
% end
% Survival_Rate=zeros(TimeLength*24+1,1);
% for t11=0:1:TimeLength*24
%     Survival_Rate(t11+1)=person_0+1-sum(C(:,t11+1)>threshold);
% end
% t=0:delt_t:delt_t*TimeLength;
% 
% 
% figure,
% for person=1:person_0+1
%     hold on;
% plot(t,C(person,:));
% end
% legend
% 
% figure,
% plot(t,Survival_Rate);
% xlabel('Time (Days)','FontWeight','Bold','FontSize',18);
% ylabel('Survival Percentage','FontWeight','Bold','FontSize',18);
% set(gca,'FontSize',14);     
% figure_FontSize=18;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',14),'FontSize',figure_FontSize);


  A1=zeros(TimeLength+1,L2+1);
  A2=zeros(TimeLength+1,L2+1);
  A3=zeros(TimeLength+1,L2+1);
  A4=zeros(TimeLength+1,L2+1);
  A5=zeros(TimeLength+1,L2+1);
  A6=zeros(TimeLength+1,L2+1);
  A7=zeros(TimeLength+1,L2+1);
  A8=zeros(TimeLength+1,L2+1);
  A9=zeros(TimeLength+1,L2+1);
  A10=zeros(TimeLength+1,L2+1);
  A11=zeros(TimeLength+1,L2+1);
  A12=zeros(TimeLength+1,L2+1);
%   A21=zeros(TimeLength+1,L2+1);
%   A31=zeros(TimeLength+1,L2+1);
%   A41=zeros(TimeLength+1,L2+1);
%   A51=zeros(TimeLength+1,L2+1);
%   A61=zeros(TimeLength+1,L2+1);

%   
  x00=0:1:L2;
  t00=0:1:TimeLength;

%   TimeLength=200;
% %   
%   L=L*x_0
global b1 b2
b1=1;b2=0;
% d_max=1;
m0=1;
  for i= 0:1:TimeLength
      for j=0:1:L2
          i1=i*delt_t;
          j1=L/L2*j;
          if j1==0
              j1=0.0000001;
          end
%           [drugt1,drugt2] = picture(j1,i1,1);
%           drugt3=Drug_34(j1,i1,1);
          drugt11=Drug2(j1,i1,1,1);
          IL411=A_Drugsimulation2(j1,i1,1,1);
          drugt13=Drug2(j1,i1,3,1);
          IL413=A_Drugsimulation2(j1,i1,3,1);
%           drug4=Drug_ceshi_jietihanshu(j1,i1);
%           drug5=Drug_5(j1*x_0,i1);
%           IL41=IL4_1(j1,i1,1);
%           IL42=IL4_2(j1,i1);
%           IL43=IL4_3(j1,i1,1);
%           IL44=A_4(j1,i1);
%           IL45=A_5(j1,i1);
%           A1(i+1,j+1)=drugt1;
%           A2(i+1,j+1)=drugt2;
%           A3(i+1,j+1)=drugt3;
%           A4(i+1,j+1)=IL41;
%           A5(i+1,j+1)=IL42;
%           A6(i+1,j+1)=IL43;
%           A7(i+1,j+1)=drug4;
%           A8(i+1,j+1)=IL44;
          A9(i+1,j+1)=drugt11;
          A10(i+1,j+1)=IL411;
          A11(i+1,j+1)=drugt13;
          A12(i+1,j+1)=IL413;
      end
  end
  
%  save('A1A3A4A6的值来画图_3_zhouqi3_0_1.mat','A1','A3','A4','A6')
%  load('A1A3A4A6的值来画图_2.mat','A1','A3','A4','A6')
   figure
%%%%% surf(x,t,u1)
pcolor(x00,t00,A1)
% title('Drug Concentration','FontWeight','Bold','FontSize',18);
title('Drug Concentration');
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:20:L2]); 
   set(gca,'xticklabel',0:20:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position');ylabel('Time (Days)');
%    set(gca,'FontWeight','Bold','FontSize',14);
   set(gca,'FontSize',16);
   
     figure
%%%%% surf(x,t,u1)
pcolor(x00,t00,A3)
title('Drug3 Concentration','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:20:L2]); 
   set(gca,'xticklabel',0:20:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',14);
     figure
%%%%% surf(x,t,u1)
pcolor(x00,t00,A4)
title('A1 Concentration','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:20:L2]); 
   set(gca,'xticklabel',0:20:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',14);
   figure
%%%%% surf(x,t,u1)
pcolor(x00,t00,A6)
title('A3 Concentration','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:20:L2]); 
   set(gca,'xticklabel',0:20:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',14);
      figure
%%%%% surf(x,t,u1)
pcolor(x00,t00,A9)
title('Drug11 Concentration','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:L2/5:L2]); 
   set(gca,'xticklabel',0:L2/5:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',14);
         figure
%%%%% surf(x,t,u1)
pcolor(x00,t00,A10)
title('A11 Concentration','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:L2/5:L2]); 
   set(gca,'xticklabel',0:L2/5:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',14);
   figure
%%%%% surf(x,t,u1)
pcolor(x00,t00,A11)
title('Drug13 Concentration','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:L2/5:L2]); 
   set(gca,'xticklabel',0:L2/5:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',14);
         figure
%%%% surf(x,t,u1)
pcolor(x00,t00,A12)
title('A13 Concentration','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0,3],'ytick',(0:0.5:3));
colorbar; 
shading interp;
   set(gca,'xtick',[0:L2/5:L2]); 
   set(gca,'xticklabel',0:L2/5:L2);
   set(gca,'ytick',0:TimeLength/5:TimeLength);  
   set(gca,'yticklabel',0:TimeLength/5:TimeLength)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',14);
  
  
  
  
%   D_d1=D_d/x_0^2;%药物的扩散率
%   L=L/x_0;
%   for i= 0:1:TimeLength
%       for j=1:1:101
%           i1=i*delt_t;
%           j1=j;
%           [drugt1,drugt2] = picture(j1,i1);
%           drugt5=Drug_34(j1,i1);
%           IL43=IL4_3(j1,i1);
%           IL41=IL4_1(j1,i1);
%           IL42=IL4_2(j1,i1);
%           A11(i+1,j)=drugt1;
%           A21(i+1,j)=drugt2;
%           A31(i+1,j)=drugt5;
%           A41(i+1,j)=IL43;
%           A51(i+1,j)=IL41;
%           A61(i+1,j)=IL42;
%       end
%   end
% 
% 
% 
% 
%
%    
%    figure
% %%%%% surf(x,t,u1)
% pcolor(x00,t00,A1)
% title('Drug concentration','FontWeight','Bold','FontSize',18);
% % colorbar('ylim',[0,3],'ytick',(0:0.5:3));
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:TimeLength/5:TimeLength);  
%    set(gca,'yticklabel',0:TimeLength/5:TimeLength)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',14);
% % 
%    figure
% % % % surf(x,t,u1)
% pcolor(x00,t00,A2)
% title('Drug 2','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength); 
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % % surf(x,t,u1)
% pcolor(x00,t00,A3)
% title('Drug 3','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength); 
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % surf(x,t,u1)
% pcolor(x00,t00,A4)
% title('IL 3','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength);  
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % surf(x,t,u1)
% pcolor(x00,t00,A5)
% title('IL 1','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength);  
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% % 
% figure
% % % surf(x,t,u1)
% pcolor(x00,t00,A6)
% title('IL 2','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength);  
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
%  figure
% %%%%% surf(x,t,u1)
% pcolor(x00,t00,A11)
% title('Drug concentration','FontWeight','Bold','FontSize',18);
% % colorbar('ylim',[0,3],'ytick',(0:0.5:3));
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:TimeLength/5:TimeLength);  
%    set(gca,'yticklabel',0:TimeLength/5:TimeLength)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',14);
% % 
%    figure
% % % % surf(x,t,u1)
% pcolor(x00,t00,A21)
% title('Drug 2','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength); 
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % % surf(x,t,u1)
% pcolor(x00,t00,A31)
% title('Drug 3','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength); 
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % surf(x,t,u1)
% pcolor(x00,t00,A41)
% title('IL 3','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength);  
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % surf(x,t,u1)
% pcolor(x00,t00,A51)
% title('IL 1','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength);  
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% % 
% figure
% % % surf(x,t,u1)
% pcolor(x00,t00,A61)
% title('IL 2','FontWeight','Bold','FontSize',18);
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:20:TimeLength);  
%    set(gca,'yticklabel',0:20:TimeLength);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
%  
%  [qq1,qq2]=meshgrid(0:10:100)
%  plot(qq1,qq2,'k',qq2,qq1,'k');
%  axis equal
%  axis([0 100 0 100])
%    set(gca,'xtick',[0:10:100]); 
%    set(gca,'xticklabel',0:10:100);
%    set(gca,'ytick',0:10:100); 
%    set(gca,'yticklabel',0:20:200)

figure
% surf(x,t,u1)
pcolor(x,t,u1)
% title('Tumor')
% title('Cancer Cell Density','FontWeight','Bold','FontSize',18);
title('Cancer Cell Density');
% colorbar('ylim',[0.3,0.7],'ytick',(0.3:0.1:0.7));
% colorbar('Ticks',0:0.1:1);
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100);
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18;
   xlabel('Position');ylabel('Time (Days)');
%    set(gca,'FontWeight','Bold','FontSize',14);
   set(gca,'FontSize',16);


figure
% surf(x,t,u2)
pcolor(x,t,u2)
title('M1 Macrophage Density');
% colorbar('ylim',[0.3,1],'ytick',(0.3:0.1:1));
colorbar; 
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L);set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);


figure
% surf(x,t,u3)
pcolor(x,t,u3)
title('M2 Macrophage Density');
colorbar;
% colorbar('ylim',[0.3,1],'ytick',(0.3:0.1:1));
 shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L);set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);


figure
% surf(x,t,u4)
pcolor(x,t,u4)
title('CSF1 Concentration');
% colorbar('ylim',[0.1,0.4],'ytick',(0.1:0.05:0.4)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L);set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);


 figure
 % surf(x,t,u5)
 pcolor(x,t,u5)
 title('EGF Concentration');
% colorbar('ylim',[0.3,0.6],'ytick',(0.3:0.05:0.6)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);
    
 figure
 % surf(x,t,u5)
 pcolor(x,t,u6)
 title('IGF1 Concentration');
% colorbar('ylim',[0,0.9],'ytick',(0:0.1:0.9));
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L);set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);
    
figure
 % surf(x,t,u5)
 pcolor(x,t,u7)
 title('EGFR Concentration');
% colorbar('ylim',[0.7,0.78],'ytick',(0.7:0.01:0.78)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L);set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);
    
figure
 % surf(x,t,u5)
 pcolor(x,t,u8)
 title('ERK Concentration');
% colorbar('ylim',[0.55,0.75],'ytick',(0.55:0.05:0.75));
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);
    
figure
 % surf(x,t,u5)
 pcolor(x,t,u9)
 title('IGF1R Concentration');
% colorbar('ylim',[0,0.025],'ytick',(0:0.005:0.025)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L);set(gca,'xticklabel',0:20:100);
   xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);
    
figure
 % surf(x,t,u5)
 pcolor(x,t,u10)
 title('AKT Concentration');
% colorbar('ylim',[0,0.45],'ytick',(0:0.05:0.45)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength);
   set(gca,'xtick',0:L/5:L);set(gca,'xticklabel',0:20:100);
    xlabel('Position');ylabel('Time (Days)');
   set(gca,'FontSize',16);
%  load('yongyaoduibi1_CSF1RI_0.5andIGF1RI_1.mat') 
%  TimeLength=1000;

   figure
   plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);
   title('Average Cancer Cell Density','FontWeight','Bold','FontSize',18)
   set(gca,'xtick',(0:delt_t*TimeLength/5:delt_t*TimeLength));
   set(gca,'xticklabel',0:TimeLength/5:TimeLength);
   set(gca,'ytick',(0:0.1:0.8));
   set(gca,'yticklabel',0:0.1:0.8)
   set(gca,'xlim',[0 delt_t*TimeLength]);
%    set(gca,'ylim',[0.25 0.7]);
   xlabel('Days','FontWeight','Bold','FontSize',18);ylabel('Average Cancer Cell Density','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',16);    

%  load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据\grid_search_CSF1RIIN_IGF1RCON_4_Jian_zhouqi3_2.mat','ObjV','parameter1','C_T_mean_1');  
%  load('D:\Matlab代码存放\2018胶质瘤\总程序\grid_search_CSF1RICON_IGF1RICON_3_lian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');    
%  load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据\grid_search_CSF1RICON_IGF1RCON_3_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');      
%  
%  load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\500天\new_grid_search_CSF1RICON_0.08_0.02_0.14_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
%  C_T_mean_11=C_T_mean_1;
% %  load('D:\Matlab代码存放\2018胶质瘤\总程序\全新的CSF1RI\grid_search_CSF1RIIN_0.01_0.01_0.1_zhouqi5.mat','ObjV','parameter1','C_T_mean_1');
% %  C_T_mean_13=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\500天\new_20_grid_search_CSF1RIIN_0.08_0.02_0.14_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
% C_T_mean_12=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\500天\new_50_grid_search_CSF1RIIN_0.08_0.02_0.14_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
% C_T_mean_13=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\500天\new_80_grid_search_CSF1RIIN_0.08_0.02_0.14_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
% C_T_mean_14=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\grid_search_CSF1RICON_0.045_0.005_0.075_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
% C_T_mean_15=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\全新的CSF1RI_EGFRI\10_10_100_grid_search_CSF1RIIN_0.09_0.01_0.15_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
% C_T_mean_16=C_T_mean_1;
% % ObjV3=ObjV;
% % ObjV2(4:5,:)=ObjV;
% % ObjV=ObjV2
% % C_T_mean_1=C_T_mean_12
% % parameter1(1,1)=0.04
% % parameter1(1,2)=0.04
% % parameter1(1,3)=0.2
% % parameter1(2,1)=0.01
% % parameter1(2,2)=0.01
% % parameter1(2,3)=0.02
% % parameter1(3,1)=0.12
% % parameter1(3,2)=0.02
% % parameter1(3,3)=0.2
% % ObjV2=zeros(15,1);
% % C_T_mean_1=C_T_mean_111;
% % C_T_mean_111(:,10:15)=C_T_mean_12;
% % ObjV2(10:15)=ObjV;
% % ObjV2=ObjV;
% % ObjV(1:9,1)=ObjV2(1:9,1)
% 
% % load('D:\Matlab代码存放\2018胶质瘤\总程序\全新的CSF1RI_EGFRI\周期为7\7_grid_search_CSF1RIIN_0.16_0.04_0.2_EGFRICON_0.01_0.01_0.02_zhouqiwei7_Jian.mat','ObjV','parameter1','C_T_mean_1');
% % C_T_mean_15=C_T_mean_1;
% 
%   load('D:\Matlab代码存放\2018胶质瘤\总程序\2月12日_1CM\new_3_grid_search_CSF1RIIN_1_EGFRI_1.mat','ObjV','parameter1','C_T_mean_1');
%  C_T_mean_11=C_T_mean_1;
%   load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\用药CSF1RI一段时间停\CSF1RI_0.1_50_EGFRI_0.01.mat','u1','C_T_mean');
%  C_T_mean_12=C_T_mean;
%  load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\用药CSF1RI一段时间停\CSF1RI_0.1_80_EGFRI_0.01.mat','u1','C_T_mean');
%  C_T_mean_13=C_T_mean;
%  load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\用药CSF1RI一段时间停\CSF1RI_0.1_100_EGFRI_0.01.mat','u1','C_T_mean');
%  C_T_mean_14=C_T_mean;
% %   load('D:\Matlab代码存放\2018胶质瘤\总程序\全新的CSF1RI_EGFRI\周期为10\grid_search_CSF1RIIN_0.5_0.5_1_EGFRICON_0.5_0.5_1_Jian_zhouqiwei10_11月16号.mat','ObjV','parameter1','C_T_mean_1');
% %  C_T_mean_16=C_T_mean_1;
% 
%   C_T_mean_55=C_T_mean_12(:,1:2,2)- C_T_mean_14
%   [aaaa1 bbb1 cccc1]=find(C_T_mean_55<0)
% % C_T_mean_13=zeros(301,3,2);
% % C_T_mean_13(:,1:2,:)=C_T_mean_1;
% %  load('D:\Matlab代码存放\2018胶质瘤\总程序\grid_search_CSF1RICON_3_lian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% %  C_T_mean_13=C_T_mean_1;
% %  load('D:\Matlab代码存放\2018胶质瘤\总程序\grid_search_CSF1RIIN_EGFRICON_4_Jian_zhouqi3_2.mat','ObjV','parameter1','C_T_mean_1');
% % C_T_mean_13(:,3,:)=C_T_mean_1;
% figure
% for i=1:1:1
% %         figure,
% %         load('wufu\grid_search_CSF1RICON_2_lian.mat','ObjV','parameter1','C_T_mean_1');
%         aa=2;bb=2;
%         figure,
% %         C_T_mean=C_T_mean_11(:,2,1);
% % C_T_mean=C_T_mean_11(:,:,i);
% C_T_mean=C_T_mean_11;
%         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% %         figure
% %         C_T_mean=C_T_mean_11(:,6);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% %         C_T_mean=C_T_mean_12(:,2,1);
% % C_T_mean=C_T_mean_12;
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% % %         figure
% % %         C_T_mean=C_T_mean_12(:,i,1);
% % %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% % % end
% % %         C_T_mean=C_T_mean_13(:,2,1);
% % C_T_mean=C_T_mean_13;
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% % %         figure
% % %         C_T_mean=C_T_mean_14(:,2,1);
% % C_T_mean=C_T_mean_14;
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% % figure
% %         C_T_mean=C_T_mean_15(:,i,1);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% % end
% %         C_T_mean=C_T_mean_13(:,12,1);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% % figure
% %         C_T_mean=C_T_mean_16(:,7,1,i);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% % end
% %         
% %         load('wufu\grid_search_CSF1RICON_EGFRICON_2_lian.mat','ObjV','parameter1','C_T_mean_1');    
% %         for ii=1:1:4
% %         C_T_mean=C_T_mean_1(:,i,ii);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% %         end
% %         C_T_mean=C_T_mean_13(:,aa,bb);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% %         C_T_mean=C_T_mean_14(:,aa,bb);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% %         C_T_mean=C_T_mean_15(:,aa,bb);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% %         C_T_mean=C_T_mean_16(:,aa,bb);
% %         plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);hold on;
% 
%         title('CSF1RI  & IGF1RI','FontWeight','Bold','FontSize',18)
% %                 title('CSF1RI','FontWeight','Bold','FontSize',18)
% %         legend('Interval period 10','Interval period 20','Interval period 30','Interval period 40','Interval period 50', ...
% %             'Interval period 60','Interval period 70','Interval period 80','Interval period 90','Interval period 100')
%         set(gca,'xtick',(0:delt_t*TimeLength/5:delt_t*TimeLength));
%         set(gca,'xticklabel',0:TimeLength/5:TimeLength);
% %         set(gca,'ytick',(0:0.1:0.9));
% %         set(gca,'yticklabel',0:0.1:0.9);
%         set(gca,'xlim',[0 delt_t*TimeLength]);
% %         set(gca,'ylim',[0 0.9]);
%         xlabel('Days','FontWeight','Bold','FontSize',18);ylabel('Average Cancer Cell Density','FontWeight','Bold','FontSize',18);
%         set(gca,'FontWeight','Bold','FontSize',16);
% %         legend('Continuous  CSF1RI Dose 1','Interval Period 3 CSF1RI Dose 1 & Continuous EGFRI Dose 1',...
% %             'Continuous  CSF1RI Dose 1 & Continuous EGFRI Dose 1','Interval Period 5 CSF1RI Dose 1 & Continuous EGFRI Dose 1',...
% %              'Interval Period 7 CSF1RI Dose 1 & Continuous EGFRI Dose 1','Interval Period 10 CSF1RI Dose 1 & Continuous EGFRI Dose 1')  
% %         legend('Continuous CSF1RI Dose 0.05','Interval Period 5 CSF1RI Dose 0.1','Continuous CSF1RI Dose 0.1')
%             
% end
% legend('20 Days CSF1RI 0.1 & Continuous EGFRI 0.01','50 Days CSF1RI 0.1 & Continuous EGFRI 0.01' ,...
%     '80 Days CSF1RI 0.1 & Continuous EGFRI 0.01','100 Days CSF1RI 0.1 & Continuous EGFRI 0.01')
% legend('Continuous CSF1RI 0.1 & Continuous IGF1RI 0.01','Interval period 20 CSF1RI 0.1 & Continuous EGFRI 0.01' ,...
%     'Interval period 50 CSF1RI 0.1 & Continuous EGFRI 0.01','Interval period 80 CSF1RI 0.1 & Continuous EGFRI 0.01')
% legend('Continuous CSF1RI 0.1 & Continuous IGF1RI 0.01','Interval period 20 CSF1RI 0.1 & Continuous IGF1RI 0.01' ,...
%    'Interval period 50 CSF1RI 0.1 & Continuous IGF1RI 0.01' ,'Interval period 80 CSF1RI 0.1 & Continuous IGF1RI 0.01')
% % , ...'Interval period 100 CSF1RI 0.5 & Continuous IGF1RI 0.01')
% legend('Interval period 100 CSF1RI 0.1 & Continuous IGF1RI 0.01','Interval period 100 CSF1RI 0.2 & Continuous IGF1RI 0.01' ,...
%     'Interval period 100 CSF1RI 0.3 & Continuous IGF1RI 0.01','Interval period 100 CSF1RI 0.4 & Continuous IGF1RI 0.01', ...
%     'Interval period 100 CSF1RI 0.5 & Continuous IGF1RI 0.01')
% legend('1','2','3','4','5' ,'6','7','8','9' ,'10','11','12','13','14','15','16','17','18','19','20')
% legend('Interval period 7 CSF1RI 0.04 & Continuous IGFRI 0.02','Interval period 7 CSF1RI 0.08 & Continuous IGFRI 0.02',...
%     'Interval period 7 CSF1RI 0.12 & Continuous IGFRI 0.02','Interval period 7 CSF1RI 0.16 & Continuous IGFRI 0.02',...
%     'Interval period 7 CSF1RI 0.2 & Continuous IGFRI 0.02')
% 
% legend('Continuous CSF1RI 0.06','Interval period 5 CSF1RI 0.12',...
%     'Continuous CSF1RI 0.12','Continuous CSF1RI 0.06 & Continuous IGFRI0.01',...
%     'Interval period 5 CSF1RI 0.12 & Continuous IGFRI 0.01','Interval period 7 CSF1RI 0.12 & Continuous IGFRI 0.01',...
%     'Continuous CSF1RI 0.12 & Continuous IGFRI 0.01')
% legend('CSF1RI 0.01->EGFRI 0.02','CSF1RI 0.02->EGFRI 0.02',...
%     'CSF1RI 0.03->EGFRI 0.02','CSF1RI 0.04->EGFRI 0.02',...
%     'CSF1RI 0.05->EGFRI 0.02','CSF1RI 0.06->EGFRI 0.02',...
%     'CSF1RI 0.07->EGFRI 0.02','CSF1RI 0.08->EGFRI 0.02',...
%     'CSF1RI 0.09->EGFRI 0.02','CSF1RI 0.10->EGFRI 0.02',...
%     'CSF1RI 0.11->EGFRI 0.02','CSF1RI 0.12->EGFRI 0.02',...
%     'CSF1RI 0.13->EGFRI 0.02','CSF1RI 0.14->EGFRI 0.02',...
%     'CSF1RI 0.15->EGFRI 0.02','CSF1RI 0.16->EGFRI 0.02',...
%     'CSF1RI 0.17->EGFRI 0.02','CSF1RI 0.18->EGFRI 0.02',...
%     'CSF1RI 0.19->EGFRI 0.02','CSF1RI 0.20->EGFRI 0.02')
% legend('Continuous  CSF1RI 0.045 & Continuous IGF1RI 0.01','Continuous  CSF1RI 0.05 & Continuous IGF1RI 0.01',...
%     'Continuous  CSF1RI 0.055 & Continuous IGF1RI 0.01','Continuous  CSF1RI 0.06 & Continuous IGF1RI 0.01',...
%     'Continuous  CSF1RI 0.065 & Continuous IGF1RI 0.01','Continuous  CSF1RI 0.07 & Continuous IGF1RI 0.01',...
%     'Continuous  CSF1RI 0.075 & Continuous IGF1RI 0.01' )
% %     ,'Continuous  CSF1RI 0.8 & Continuous IGF1RI 0.01',...
% %     'Continuous  CSF1RI 0.9 & Continuous IGF1RI 0.01')
% legend('Continuous  CSF1RI 0.045 & Continuous EGFRI 0.01','Continuous  CSF1RI 0.05 & Continuous EGFRI 0.01',...
%     'Continuous  CSF1RI 0.055 & Continuous EGFRI 0.01','Continuous  CSF1RI 0.06 & Continuous EGFRI 0.01',...
%     'Continuous  CSF1RI 0.065 & Continuous EGFRI 0.01','Continuous  CSF1RI 0.07 & Continuous EGFRI 0.01',...
%     'Continuous  CSF1RI 0.075 & Continuous EGFRI 0.01')
% % ,'Continuous  CSF1RI 0.8 & Continuous EGFRI 0.01',...
% %     'Continuous  CSF1RI 0.9 & Continuous EGFRI 0.01')
% %     'Continuous  CSF1RI 0.11','Continuous  CSF1RI 0.12',...
% %     'Continuous  CSF1RI 0.13','Continuous  CSF1RI 0.14',...
% %     'Continuous  CSF1RI 0.15','Continuous  CSF1RI 0.16',...
% %     'Continuous  CSF1RI 0.17','Continuous  CSF1RI 0.18',...
% %     'Continuous  CSF1RI 0.19','Continuous  CSF1RI 0.2');
% 
     
   
end


function [c,f1,s1] = pdex4pde(x,t,u,DuDx)
global A_T_EGF A_T_IGF1 A_M1 A_M2 epsilon r_T a_ERK K1 K2 K31 K32 K4 a_AKT  d_T d_TM1 a_M2M1 a_CSF1 a_A B_1 B_2 b_EGF
global CSF1_max EGF_max IGF1_max
global ERK AKT 
global x_0
global Ptest1
global delt_t
global select_CSF1R_I dmax1
% global drug
% Parameters_nondimensional_4_gai_jin_fangcheng_4
c = [1; 1; 1; 1; 1; 1];  
% EGFR_I=0;
% IGF1R_I=0;
% x1=x/x_0;

% bbb=abs(u)<=0.005;
%     u(bbb)=0.005;
t
f1 = [1,0,0,0,-A_T_EGF*u(1),-A_T_IGF1*u(1);%A_T是新的A1了
      0,1,0,-A_M1*u(2),0,0;%A是新置换的巨噬细胞的参数
      0,0,1,-A_M2*u(3),0,0;
      0,0,0,1/epsilon,0,0;
      0,0,0,0,1/epsilon,0;
      0,0,0,0,0,1/epsilon]*DuDx; 

% f = [1,0,0,0,-A_T_EGF*u(1),-A_T_IGF1*u(1),0,0,0,0;%A_T是新的A1了
%     0,1,0,-A_M1*u(2),0,0,0,0,0,0;%A是新置换的巨噬细胞的参数
%     0,0,1,-A_M2*u(3),0,0,0,0,0,0;
%     0,0,0,1/epsilon,0,0,0,0,0,0;
%     0,0,0,0,1/epsilon,0,0,0,0,0;
%     0,0,0,0,0,1/epsilon,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0]*DuDx; 

% f = [0,0,0,0;
%     0,0,0,0;
%     0,0,0,0;
%     0,0,0,0]*DuDx; 
% % 
%      drug=0;%边界条件药物为0
%      IL4=0; 
% if t>=delt_t*25
%      t1=t-delt_t*25;
%      [drug,buyao] = picture(x,t1); %边界条件一药物的值    
%      IL4=IL4_1(x,t1);
% else
%     drug=0;
%     IL4=0; 
% end
%      [drug,buyao] = picture(x,t); %边界条件一药物的值    
%      IL4=IL4_1(x,t);
%      [buyao,drug] = picture(x,t); %边界条件二的时候药物的值
%      IL4=IL4_2(x,t);
%      drug=Drug_34(x,t);%边界条件三的时候药物的值 
%      IL4=IL4_3(x,t);

% %调试
% %      R1=r_T*(1+a_r*u(5)./(K1+u(5))+a_r*u(6)./(K4+u(6)));  % growth rate
%      R1=r_T*(1+a_ERK*value_ERK/(K5+value_ERK)+a_AKT*value_AKT/(K6+value_AKT));
%      D1=d_T*(1+d_TM1*u(2));  % death rate    
% %      R3=a_M2M1*u(3)-u(2).*(a_CSF1*u(4)./(K21+K22*drug+u(4))+a_A*IL4/(K3+IL4));
%      R3=a_M2M1*u(3)-u(2)*(a_CSF1*u(4)/(K21+K22*drug+u(4))+a_A*IL4/(K3+IL4));
%      bs=u(3)*IL4-u(6);
%     if bs<0
%         bs=0;
%     end
     drug_CSF1R_I=Drug2(x,t,select_CSF1R_I,dmax1,0);
     IL4=A_Drugsimulation2(x,t,select_CSF1R_I,dmax1,0);
%调试
%      R1=r_T*(1+a_r*u(5)./(K1+u(5))+a_r*u(6)./(K4+u(6)));  % growth rate
     R1=r_T*(1+a_ERK*u(8)/(K1+u(8))+a_AKT*u(10)/(K2+u(10)));
     D1=d_T*(1+d_TM1*u(2));  % death rate    
%      R3=a_M2M1*u(3)-u(2).*(a_CSF1*u(4)./(K21+K22*drug_CSF1R_I+u(4))+a_A*IL4/(K3+IL4));
     R3=a_M2M1*u(3)-u(2)*(a_CSF1*u(4)/(K31+K32*drug_CSF1R_I+u(4))+a_A*IL4/(K4+IL4));
%      bs=u(3)*IL4-u(6);
%     if bs<0
%         bs=0;
%     end

%     R7=V1*EGF/(K11+EGF)/(1+u(2)/K12)/(1+EGFR_I/K13)*(1-u(1))-d1*u(1);
%     R8=a*(1+V21*u(1)^n/(K21^n+u(1)^n))*(1+V22*u(3)/(K22+u(3)))/(1+u(4)/K23)*(1-u(2))-d2*u(2);
%     R9=V3*IGF1/(K31+IGF1)/(1+u(2)/K32)/(1+IGF1R_I/K33)*(1-u(3))-d3*u(3);
%     R10=V4*u(1)/(K41+u(1))*u(3)/(K42+u(3))*(1-u(4))-d4*u(4);
% 新的
%  s = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)-u(4));1/epsilon*B_2*(u(3).*u(4)./(K21+K22*drug_CSF1R_I+u(4))-u(5));1/epsilon*bs;
%      R7;R8;R9;R10]; 
%  s1 = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)*(1-u(4)/CSF1_max)-u(4));1/epsilon*B_2*(b_EGF*u(3).*u(4)./(K210+K220*drug_CSF1R_I+u(4))*(1-u(5)/EGF_max)-u(5));1/epsilon*(u(3)*IL4*(1-u(6)/IGF1_max)-u(6))];
s1 = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)*(1-u(4)/CSF1_max)-u(4));1/epsilon*B_2*(b_EGF*u(3)*(1-u(5)/EGF_max)-u(5));1/epsilon*(u(3)*IL4*(1-u(6)/IGF1_max)-u(6))]; 
%  s2 = [R7;R8;R9;R10]; 
end
% --------------------------------------------------------------
function u0 = pdex4ic(x)
global L
global C_T_0 C_M1_0 C_M2_0 CSF1_0 EGF_0 IGF1_0
global x_0
x1=x/x_0;
% a=floor(x/x_0);
% b=x-x_0*a;
% a2=a+1;
% if b~=0
%     value1=C_T_0(1,a2)+b/x_0*(C_T_0(1,a2+1)-C_T_0(1,a2));
%     value2=C_M1_0(1,a2)+b/x_0*(C_M1_0(1,a2+1)-C_M1_0(1,a2));
%     value3=C_M2_0(1,a2)+b/x_0*(C_M2_0(1,a2+1)-C_M2_0(1,a2));
%     value4=CSF1_0(1,a2)+b/x_0*(CSF1_0(1,a2+1)-CSF1_0(1,a2));
%     value5=EGF_0(1,a2)+b/x_0*(EGF_0(1,a2+1)-EGF_0(1,a2));
%     value6=IGF1_0(1,a2)+b/x_0*(IGF1_0(1,a2+1)-IGF1_0(1,a2));
% else
%     value1=C_T_0(1,a2);
%     value2=C_M1_0(1,a2);
%     value3=C_M2_0(1,a2);
%     value4=CSF1_0(1,a2);
%     value5=EGF_0(1,a2);
%     value6=IGF1_0(1,a2);
% end
ExpData;
     EGFR_0=pY1068_EGFR(1,1);%initial value
     ERK_0=p_ERK(1,1);
     IGF1R_0=p_IGF1Rbeta(1,1);
     AKT_0=pS473_AKT(1,1);
% u0 = [0.5.*Heaviside(11-x).*(1+0.1.*randn(1,length(x))); 100/100.*(1+0.1.*randn(1,length(x)))+0.0; 1.*(1+0.01.*randn(1,length(x))); 1.*(1+0.01.*randn(1,length(x))); d_max*x/100];    

% u0 = [(1-x/100)*(0.1+0*0.1.*rand(1,length(x)))+0.5;x/100*(0.05+0*0.1.*randn(1,length(x)))+0.12; x/100.*(0.05+0*0.1.*randn(1,length(x)))+0.88; (1-x/100)*(0.1+0*0.1.*rand(1,length(x)))+0.7; (1-x/100)*(0.2+0*0.1.*rand(1,length(x)))+0.3;0;
%     0.7;0.72;0.01;0.01];   %用药的时候
u0 = [(1-x/L)*(0.1+0*0.1.*rand(1,length(L)))+0.5;x/L*(0.05+0*0.1.*randn(1,length(L)))+0.12; x/L.*(0.05+0*0.1.*randn(1,length(x)))+0.88; (1-x/L)*(0.1+0*0.1.*rand(1,length(x)))+0.4; (1-x/L)*(0.1+0*0.1.*rand(1,length(x)))+0.4;0.01;
    0.7;0.72;0.01;0.01];   %用药的时候   改了CSF1

% u0 = [(1-x/100)*(0.1)+0.08;x/100*(0.05)+0.12; x/100.*(0.05)+0.88; (1-x/100)*(0.1+0*0.1.*rand(1,length(x)))+0.7; (1-x/100)*(0.2+0*0.1.*rand(1,length(x)))+0.3;0;
%     0.7;0.72;0.01;0.01];   %从小开始生长



end


% --------------------------------------------------------------

function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
pl = [0; 0; 0; 0; 0; 0];                               
ql = [1; 1; 1; 1; 1; 1];
pr = [0; 0; 0; 0; 0; 0];
%pr = [0; 0; 0; 0; ur(5)-dnt];
%pr = [0; 0; 0; 0; ur(5)-d_max/2*(sin((t/delt_t)*pi/100)+2)];
%pr = [0; 0; 0; 0; ur(5)-(3-rem(floor(t/(delt_t*30)),3))];%忘记了是哪种情况了
qr = [1; 1; 1; 1; 1; 1];
% pl = [0; 0; 0; 0];                               
% ql = [0; 0; 0; 0];       
% pr = [0; 0; 0; 0];       
% %pr = [0; 0; 0; 0; ur(5)-dnt];
% %pr = [0; 0; 0; 0; ur(5)-d_max/2*(sin((t/delt_t)*pi/100)+2)];
% %pr = [0; 0; 0; 0; ur(5)-(3-rem(floor(t/(delt_t*30)),3))];%忘记了是哪种情况了
% qr = [0; 0; 0; 0];       
end

function [drugt1,drugt2] = picture(x,t,dmax)
global d_max
global eta1
global D_d1
global L   %x的最大值
global delt_t
global TimeLength
d_initial=d_max/L;
w=pi/(TimeLength*delt_t/20);
   d11=dmax;
   d110=dmax;
   d22=d_max/2*(sin(w*t)+2);
drugt1=0;
drugt2=0;
drugt11=0;
drugt21=0;
for n=1:1:100
    Tn1=((2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+  d110*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t)))*sin(n*pi*x/L);
    drugt11=drugt11+Tn1; 
    Tn2=(((2*L^2*(-1)^(n+1)/(n*pi)+4*L^2*(-1)^n/(n*pi)^3-4*L^2/(n*pi)^3)*d_initial+2*L*(-1)^n/(n*pi)*d_max)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t))/(n*pi*D_d1*(n*pi/L)^2)+d_max*w*2*L*(-1)^n*(w*sin(w*t)+D_d1*(n*pi/L)^2*(cos(w*t)-exp(-D_d1*(n*pi/L)^2*t)))/(2*n*pi*(w^2+D_d1^2*(n*pi/L)^4)))*sin(n*pi*x/L);
    drugt21=drugt21+Tn2;
end
drugt1=1/x*drugt11+d11;
drugt2=1/x*drugt21+d22;
end



%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用的sin函数来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分


  
%下面是计算在位置x和时间t时由药物CSF1R_I产生的药物产生的IL4的浓度

%%药物边界条件为常量3的时候产生的IL4的量
%第一部分
function q = TnIL4_1(x,t,n,dmax)
% global d_max
global D_d1;
global eta1;
global L
d_initial=dmax/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
% q2=0;%因为又边界输入为常数3，导数值为0；
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+dmax*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/(D_d1*(n*pi/L)^2)+eta1*2*(-1)^n*t/(D_d1*(n*pi/L)^3)+eta1*2*(-1)^n*(exp(-D_d1*(n*pi/L)^2*t)-1)/(D_d1^2*(n*pi/L)^5);

end
%第二部分
function q = IL4_1(x,t,dmax)
% global d_max
global A_0
global L
global B_3  %IL4的分泌速率 
q0=0;
for n=1:1:100
    q0=q0+sin(n*pi/L*x)*TnIL4_1(x,t,n,dmax);
end
q=A_0+B_3/x*q0+B_3*dmax*t;
end

%%药物边界条件为d22=d_max/2*(sin(w*t)+2)=3/2*(sin(w*t)+2)的时候产生的IL4的量
%第一部分
function q5 = TnIL4_2(t,n,dmax)
% global d_max
global D_d1;
global eta1;
global L
global delt_t
global TimeLength
w=pi/(TimeLength*delt_t/20);
d_initial=dmax/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
% q2=0;%因为又边界输入为常数3，导数值为0；
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+dmax*(2*L)/n/pi*(-1)^n;
q1=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/(D_d1*(n*pi/L)^2);
q2=eta1*2*(-1)^n*t/(D_d1*(n*pi/L)^3);
q3=eta1*2*(-1)^n/(D_d1^2*(n*pi/L)^5)*(exp(-D_d1*(n*pi/L)^2*t)-1);
q4=dmax*w*L/n/pi*(-1)^n*(D_d1*(n*pi/L)^2*sin(w*t)/w+exp(-D_d1*(n*pi/L)^2*t)-cos(w*t))/(w^2+D_d1^2*(n*pi/L)^4);
q5=q1+q2+q3+q4;
end
%第二部分
function q = IL4_2(x,t,dmax)
% global d_max
global A_0
global L
global B_3  %IL4的分泌速率 
global TimeLength
global delt_t
q0=0;
w=pi/(TimeLength*delt_t/20);
for n=1:1:100
    q0=q0+sin(n*pi*x/L)*TnIL4_2(t,n);
end
q=A_0+B_3/x*q0+B_3*dmax/2/w*(1-cos(w*t))+B_3*dmax*t;
end

%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分
%积分函数一，对应h1(t,t1)  t1是to
function q4=h_1_ceshiair1(k,t1,t,n,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*exp(-D_d1*(n*pi/L)^2*(t-t1))*(2*m0/d^3*(t1-(k+1/2)*T)^2+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d))*(t1-(k+1/2)*T));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*exp(-D_d1*(n*pi/L)^2*(t-t1))*(8*m0/d^3*(t1-(k+1/2)*T)+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d)));
q3=24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(t-t1));
q4=q1+q2+q3;
end
%积分函数二，对应h1(t,t1)
function q4=h_2_ceshiair2(k,t1,t,n,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*exp(-D_d1*(n*pi/L)^2*(t-t1))*(-2*m0/d^3*(t1-(k+1)*T+d)^2+2*m0/d^2*(1-2/d*(t1-(k+1)*T))*(t1-(k+1)*T+d));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*exp(-D_d1*(n*pi/L)^2*(t-t1))*(-8*m0/d^3*(t1-(k+1)*T+d)+2*m0/d^2*(1-2/d*(t1-(k+1)*T)));
q3=-24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(t-t1));
q4=q1+q2+q3;
end
%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q = Tn34(t,n,m0)
global D_d1;
global eta1;
global L
global period1
global inter1
% global m0
T=2*period1;
d=inter1;
q=0;
q1=0;
q2=0;
q3=0;
q4=0;
d_0=m0;%初始条件为d_initial*r，应该要除以L
d_initial=d_0/L;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
for i=0:1:(d2-1)
    t1=(i+1/2)*T-d;
    q1=h_1_ceshiair1(i,t1,t,n,m0);
    t1=(i+1/2)*T;
    q2=h_1_ceshiair1(i,t1,t,n,m0);
    q3=q2-q1;
    q4=q4+q3;
    t1=(i+1)*T-d;
    q1=h_2_ceshiair2(i,t1,t,n,m0);
    t1=(i+1)*T;
    q2=h_2_ceshiair2(i,t1,t,n,m0);
    q3=q2-q1;
    q4=q4+q3;
end
if d3==0 && t0>=(period1-inter1)
    t1=(d2+1/2)*T-d;
    q1=h_1_ceshiair1(d2,t1,t,n,m0);
    t1=t;
    q2=h_1_ceshiair1(d2,t1,t,n,m0);
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==1 && t0<(period1-inter1)
    t1=(d2+1/2)*T-d;
    q1=h_1_ceshiair1(d2,t1,t,n,m0);
    t1=(d2+1/2)*T;
    q2=h_1_ceshiair1(d2,t1,t,n,m0);
    q3=q2-q1;
    q4=q4+q3; 
end 
if d3==1 && t0>=(period1-inter1)
    t1=(d2+1/2)*T-d;
    q1=h_1_ceshiair1(d2,t1,t,n,m0);
    t1=(d2+1/2)*T;
    q2=h_1_ceshiair1(d2,t1,t,n,m0);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1)*T-d;
    q1=h_2_ceshiair2(d2,t1,t,n,m0);
    t1=t;
    q2=h_2_ceshiair2(d2,t1,t,n,m0);
    q3=q2-q1;
    q4=q4+q3;
end 
q=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t))+q4;
end
%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q = Drug_34(x,t,m0) 
global L
global period1
global inter1
% global m0
q0=0;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
xk=(d1+1)*period1-inter1;
xk1=(d1+1)*period1;
if d3==0 && t0<=(period1-inter1)
    dnt=m0;
end
if d3==0 && t0>(period1-inter1)
    yk=m0;
    %yk1=0;
    dnt=(1+2*(t-xk)/(xk1-xk))*((t-xk1)/(xk-xk1))^2*yk+0+0+0;   
end
if d3==1 && t0<=(period1-inter1)
    dnt=0;
end
if d3==1 && t0>(period1-inter1)
    %yk=0;
    yk1=m0;
    dnt=0+(1+2*(t-xk1)/(xk-xk1))*((t-xk)/(xk1-xk))^2*yk1+0+0;
    
end
    for n=1:1:100
    q0=q0+Tn34(t,n,m0)*sin(n*pi*x/L);  
    end
q=q0/x+dnt;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数三
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%第一部分
function q6=h_3_ceshiair3(k,s,n,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*(2*m0/3/d^3*(s-(k+1/2)*T)^3+2*m0/d^2*(s^2/2+2/d*(s^3/3-((k+1/2)*T-d)/2*s^2))-2*m0*(k+1/2)*T/d^2*(s+1/d*(s-(k+1/2)*T+d)^2));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*(4*m0/d^3*(s-(k+1/2)*T)^2+2*m0/d^2*(s+1/d*(s-(k+1/2)*T+d)^2));
q3=24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*s;
q4=12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q5=24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q6=q1+q2+q3+q4+q5;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数四
function q6=h_4_ceshiair4(k,s,n,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*(-2*m0/3/d^3*(s-(k+1)*T+d)^3+2*m0/d^2*(s^2/2-2/d*(s^3/3-(k+1)*T/2*s^2))-2*m0*((k+1)*T-d)/d^2*(s-1/d*(s-(k+1)*T)^2));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*(-4*m0/d^3*(s-(k+1)*T+d)^2+2*m0/d^2*(s-1/d*(s-(k+1)*T)^2));
q3=-24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*s;
q4=-12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q5=-24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q6=q1+q2+q3+q4+q5;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数五
%对应的是对gk1函数的积分
%积分函数五是因为fn包含边界函数d'(t)是阶梯函数，分段积分完成了后是一个定的关于t的积分，再对t积分后和积分函数三不一样，然后才是IL的值的一部分
function q5=h_5_ceshiair5(k,s,n,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T));
q2=12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q3=-24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T));
q4=24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数六
%对应的是对gk2函数的积分
function q5=h_6_ceshiair6(k,s,n,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=-12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T));
q2=-12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q3=24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T));
q4=-24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数七
%对应的是对d(t)阶梯函数第一部分的积分
function q1=h_7_ceshiair7(s,m0)
% global m0
q1=m0*s;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数八
%对应的是对d(t)阶梯函数第二部分的积分
function q5=h_8_ceshiair8(k,s,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=m0/d^2/2/d*s^4;
q2=m0/d^2/3*(3-2/d*(k+1/2)*T-4/d*(k+1/2)*T)*s^3;
q3=m0/d^2/2*(-6+6/d*(k+1/2)*T)*(k+1/2)*T*s^2;
q4=m0/d^2*(3-2/d*(k+1/2)*T)*(k+1/2)^2*T^2*s;
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数八
%对应的是对d(t)阶梯函数第二部分的积分
function q5=h_9_ceshiair9(k,s,m0)
global D_d1
global L
global inter1
global period1
% global m0
T=2*period1;
d=inter1;
q1=-m0/d^2/2/d*s^4;
q2=m0/d^2/3*(1+2/d*(k+1)*T-4/d*(d-(k+1)*T))*s^3;
q3=m0/d^2/2*(-2/d*(d-(k+1)*T)^2+2*(1+2/d*(k+1)*T)*(d-(k+1)*T))*s^2;
q4=m0/d^2*(1+2/d*(k+1)*T)*(d-(k+1)*T)^2*s;
q5=q1+q2+q3+q4;
end
%%药物边界条件为阶梯函数6和0的时候产生的IL4的量，中间用埃尔米插值连接
%第二部分
function q = TnIL4_3(t,n,m0)

global D_d1;
global eta1;
global L
global period1  %period1为一个周期
global inter1
% global m0
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q=0;
q4=0;
d_0=m0;%药物右边界条件，持续输入
d_initial=d_0/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%积分函数五是因为fn包含边界函数d'(t)是阶梯函数，分段积分完成了后是一个定的关于t的积分，再对t积分后和积分函数三不一样，然后才是IL的值的一部分
for i=0:1:(d2-1)
    s=(i+1/2)*T-d;
    q1=h_3_ceshiair3(i,s,n,m0);%调用积分函数三
    s=(i+1/2)*T;
    q2=h_3_ceshiair3(i,s,n,m0);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T;
    q1=h_5_ceshiair5(i,s,n,m0);%调用积分函数五
    s=t;
    q2=h_5_ceshiair5(i,s,n,m0);%调用积分函数五
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T-d;
    q1=h_4_ceshiair4(i,s,n,m0);%调用积分函数四
    s=(i+1)*T;
    q2=h_4_ceshiair4(i,s,n,m0);%调用积分函数四
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T;
    q1=h_6_ceshiair6(i,s,n,m0);%调用积分函数六
    s=t;
    q2=h_6_ceshiair6(i,s,n,m0);%调用积分函数六
    q3=q2-q1;
    q4=q4+q3;
end
% for i=0:1:(d2-1)
%     s=(i+1/2)*T-d;
%     q1=h_3_ceshiair3(i,s,n);%调用积分函数三
%     s=(i+1/2)*T;
%     q2=h_3_ceshiair3(i,s,n);%调用积分函数三
%     q3=q2-q1;
%     q4=q4+q3;
%     s=(i+1)*T-d;
%     q1=h_4_ceshiair4(i,s,n);%调用积分函数四
%     s=(i+1)*T;
%     q2=h_4_ceshiair4(i,s,n);%调用积分函数四
%     q3=q2-q1;
%     q4=q4+q3;
% end

%判断t最后走到第几个周期积分对应着对应h1(s),h2(s)的积分
if d3==0 && t0>=(period1-inter1)
    s=(d2+1/2)*T-d;
    q1=h_3_ceshiair3(d2,s,n,m0);%调用积分函数三
    s=t;
    q2=h_3_ceshiair3(d2,s,n,m0);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==1 && t0<(period1-inter1)
    s=(d2+1/2)*T-d;
    q1=h_3_ceshiair3(d2,s,n,m0);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_3_ceshiair3(d2,s,n,m0);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T;
    q1=h_5_ceshiair5(d2,s,n,m0);%调用积分函数五
    s=t;
    q2=h_5_ceshiair5(d2,s,n,m0);%调用积分函数五
    q3=q2-q1;
    q4=q4+q3;    
end 
if d3==1 && t0>=(period1-inter1)
    s=(d2+1/2)*T-d;
    q1=h_3_ceshiair3(d2,s,n,m0);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_3_ceshiair3(d2,s,n,m0);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T;
    q1=h_5_ceshiair5(d2,s,n,m0);%调用积分函数五
    s=t;
    q2=h_5_ceshiair5(d2,s,n,m0);%调用积分函数五
    q3=q2-q1;
    q4=q4+q3;    
    s=(d2+1)*T-d;
    q1=h_4_ceshiair4(d2,s,n,m0);%调用积分函数四
    s=t;
    q2=h_4_ceshiair4(d2,s,n,m0);%调用积分函数四
    q3=q2-q1;
    q4=q4+q3;
end 

fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q5=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/D_d1/(n*pi/L)^2;
q6=eta1*2*(-1)^n/D_d1/(n*pi/L)^3*t;
q7=eta1*2*(-1)^n/D_d1^2/(n*pi/L)^5*(exp(-D_d1*(n*pi/L)^2*t)-1);
q=q4+q5+q6+q7;
end
%第三部分
function q = IL4_3(x,t,m0)
global A_0
global L
global B_3  %IL4的分泌速率 
global period1  %period1为一个周期
global inter1
% global m0
d_0=m0;%初始右边界药物的条件
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q0=0;
q1=0;
q2=0;
q3=0;
q4=0;
for k=0:1:d2-1
    s=k*T;
    q1=h_7_ceshiair7(s,m0);%调用积分函数七
    s=(k+1/2)*T-d;
    q2=h_7_ceshiair7(s,m0);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;
    s=(k+1/2)*T-d;
    q1=h_8_ceshiair8(k,s,m0);%调用积分函数八
    s=(k+1/2)*T;
    q2=h_8_ceshiair8(k,s,m0);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3;
    s=(k+1)*T-d;
    q1=h_9_ceshiair9(k,s,m0);%调用积分函数九
    s=(k+1)*T;
    q2=h_9_ceshiair9(k,s,m0);%调用积分函数九
    q3=q2-q1;
    q4=q4+q3;    
end
%判断t最后走到第几个周期积分对应着对应d(t)的最后一部分的积分
if d3==0 && t0<=(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s,m0);%调用积分函数七
    s=t;
    q2=h_7_ceshiair7(s,m0);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==0 && t0>(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s,m0);%调用积分函数七
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s,m0);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s,m0);%调用积分函数八
    s=t;
    q2=h_8_ceshiair8(d2,s,m0);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3;    
end 
if d3==1 && t0<=(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s,m0);%调用积分函数七
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s,m0);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s,m0);%调用积分函数八
    s=(d2+1/2)*T;
    q2=h_8_ceshiair8(d2,s,m0);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3; 
end 
if d3==1 && t0>(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s,m0);%调用积分函数七
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s,m0);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s,m0);%调用积分函数八
    s=(d2+1/2)*T;
    q2=h_8_ceshiair8(d2,s,m0);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1)*T-d;
    q1=h_9_ceshiair9(d2,s,m0);%调用积分函数九
    s=t;
    q2=h_9_ceshiair9(d2,s,m0);%调用积分函数九
    q3=q2-q1;
    q4=q4+q3; 
end 
for n=1:1:100
    q0=q0+TnIL4_3(t,n,m0)*sin(n*pi*x/L);
end
% q=A_0+B_3/x*q0+B_3*d_0*t+B_3*q4;
q=A_0+B_3/x*q0+B_3*q4;
end





