% function  parameter1= jianduanyongyao(a1,a2)
function  [parameter1,C_T_mean_1] = jianduanyongyao(a1,a2)

tic
clc
clear all
close all

addpath('/home/lhf/matlab code/funfun')%add operation path
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
global dmax1 dmax2 dmax3 
global select_CSF1R_I select_EGFR_I select_IGF1R_I 
Parameters_nondimensional_5_gai_jin_fangcheng_5
% load('D:\Matlab代码存放\2018胶质瘤\总程序\canshu114-24hours-ode15s.mat','p');
load('canshu1-24hours-ode15s.mat','p');
Ptest1=p(1,1:22);
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
period=40;           %小周期，例如20天用药量6，20天用药量0
period1=period*delt_t;  %真正周期
inter=1;             %间隔
inter1=inter*delt_t;   %真正间隔
eta1=eta_d;
D_d1=D_d;%药物的扩散率

if a2==1
    select_EGFR_I=1;
    select_IGF1R_I=0;
    dmax2=1;
    dmax3=0;
end
if a2==2
    select_EGFR_I=0;
    select_IGF1R_I=1;
    dmax2=0;
    dmax3=1;
end

select_CSF1R_I=3;
select_EGFR_I=1;
select_IGF1R_I=0;


% figure,
global time
time=0;
ilin=1;
iilin=5;
% iiilin=5;
ObjV=zeros(ilin,iilin);
% ObjV=zeros(ilin,1);
C_T_mean_1=zeros(TimeLength+1,ilin,iilin);
% C_T_mean_1=zeros(TimeLength+1,ilin);
for i=1:1:ilin
    for ii=1:1:iilin
%         for iii=1:1:iiilin
%             dmax1=i*0.02+0.06;
            dmax1=a1*i;
%             dmax2=1;
            




            period=ii*10;
            period1=period*delt_t;
            options=odeset('reltol',1e-3);%精度函数不能少
            % options=odeset('reltol',1);%精度函数不能少
            m = 2;%用matlab解需要设置的参数
            sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);
            if length(sol)<TimeLength+1
                u1=100*ones(TimeLength+1,101);
                C_T_mean=sum(u1(:,:),2)/101; 
                C_T_mean_end=100;   
            else
                u1 = sol(:,:,1);
                C_T_mean=sum(u1(:,:),2)/101; 
                C_T_mean_end=C_T_mean(TimeLength+1,1); 
            end
            if isempty(find(C_T_mean<0))
            else
                u1=100*ones(TimeLength+1,101);
                C_T_mean=sum(u1(:,:),2)/101; 
                C_T_mean_end=100;  
            end
            C_T_mean_1(:,i,ii)=C_T_mean;
            ObjV(i,ii)=sum(C_T_mean)+dmax1*(TimeLength+1)/2;   % y is simulated data; y0 is experimental data
%         end
    end
end
parameter1=zeros(2,3);
parameter1(1,1)=1;parameter1(1,2)=0.1;parameter1(1,3)=1;
parameter1(2,1)=1;parameter1(2,2)=1;parameter1(2,3)=1;
parameter1(3,1)=10;parameter1(3,2)=10;parameter1(3,3)=100;
% save('grid_search_CSF1RICON_EGFRICON_3_lian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1','u1');
% load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\80_90_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
filenm=['quanxin2/new_10_10_50_grid_search_CSF1RIIN_' num2str(a1) '_IGF1RICON_1.mat'];
save(filenm,'ObjV','parameter1','C_T_mean_1');


% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_0.2_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_1=ObjV;
% parameter1_1=parameter1;
% C_T_mean_1_1=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_0.4_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_2=ObjV;
% parameter1_2=parameter1;
% C_T_mean_1_2=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_0.6_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_3=ObjV;
% parameter1_3=parameter1;
% C_T_mean_1_3=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_0.8_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_4=ObjV;
% parameter1_4=parameter1;
% C_T_mean_1_4=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_1_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_5=ObjV;
% parameter1_5=parameter1;
% C_T_mean_1_5=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_1.2_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_6=ObjV;
% parameter1_6=parameter1;
% C_T_mean_1_6=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_1.4_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_7=ObjV;
% parameter1_7=parameter1;
% C_T_mean_1_7=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_1.6_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_8=ObjV;
% parameter1_8=parameter1;
% C_T_mean_1_8=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_1.8_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_9=ObjV;
% parameter1_9=parameter1;
% C_T_mean_1_9=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_2_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_10=ObjV;
% parameter1_10=parameter1;
% C_T_mean_1_10=C_T_mean_1;
% 
% 
% ObjV=zeros(10,10);
% parameter1=zeros(2,3);
% C_T_mean_1=zeros(301,10,10);
% 
% ObjV(1,:)=ObjV_1;ObjV(2,:)=ObjV_2;
% ObjV(3,:)=ObjV_3;ObjV(4,:)=ObjV_4;
% ObjV(5,:)=ObjV_5;ObjV(6,:)=ObjV_6;
% ObjV(7,:)=ObjV_7;ObjV(8,:)=ObjV_8;
% ObjV(9,:)=ObjV_9;ObjV(10,:)=ObjV_10;
% parameter1(1,1)=0.2,parameter1(1,2)=0.2,parameter1(1,3)=2;
% parameter1(2,1)=0.2,parameter1(2,2)=0.2,parameter1(2,3)=2;
% C_T_mean_1(:,1,:)=C_T_mean_1_1;C_T_mean_1(:,2,:)=C_T_mean_1_2;
% C_T_mean_1(:,3,:)=C_T_mean_1_3;C_T_mean_1(:,4,:)=C_T_mean_1_4;
% C_T_mean_1(:,5,:)=C_T_mean_1_5;C_T_mean_1(:,6,:)=C_T_mean_1_6;
% C_T_mean_1(:,7,:)=C_T_mean_1_7;C_T_mean_1(:,8,:)=C_T_mean_1_8;
% C_T_mean_1(:,9,:)=C_T_mean_1_9;C_T_mean_1(:,10,:)=C_T_mean_1_10;
% 
% C_T_mean_1_min_value=zeros(10,10);
% C_T_mean_1_min_position=zeros(10,10);
% C_T_mean_1_end_value=zeros(10,10);
% for j=1:1:10
%     for jj=1:1:10
%         [aa1,bb1]=min(C_T_mean_1(:,j,jj));
%         C_T_mean_1_min_value(j,jj)=aa1;
%         C_T_mean_1_min_position(j,jj)=bb1;
%         C_T_mean_1_end_value(j,jj)=C_T_mean_1(301,j,jj);
%     end
% end
% C_T_mean_1_reappear=zeros(10,10);
% for j=1:1:10
%     for jj=1:1:10
%         cc1=C_T_mean_1_min_position(j,jj);
%         cc2=find(C_T_mean_1(cc1:301,j,jj)>=0.6);
%         if isempty(cc2)
%             cccc=0
%         else
%             C_T_mean_1_reappear(j,jj)=C_T_mean_1_min_position(j,jj)+cc2(1)-1;
%         end
%     end
% end
% % save('grid_search_CSF1RIIN_0.2_0.2_2_IGF1RICON_0.2_0.2_2_Jian_zhouqi3_11月10号.mat','ObjV','parameter1','C_T_mean_1');
% % save('grid_search_CSF1RICON_0.5_0.5_1_lian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% % load('grid_search_CSF1RIIN_0.2_IGF1RICON_1_0.2_1.2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');

toc
end