load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\10_20_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
ObjV_1=ObjV;
parameter1_1=parameter1;
C_T_mean_1_1=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\30_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_2=ObjV;
parameter1_2=parameter1;
C_T_mean_1_2=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\40_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
ObjV_3=ObjV;
parameter1_3=parameter1;
C_T_mean_1_3=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\50_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_4=ObjV;
parameter1_4=parameter1;
C_T_mean_1_4=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\60_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
ObjV_5=ObjV;
parameter1_5=parameter1;
C_T_mean_1_5=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\70_grid_search_CSparameter1F1RIIN_0.1_0.1_2_EGFRICON_0.01_0.01_0.05_zhouqiwei70_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_6=ObjV;
parameter1_6=parameter1;
C_T_mean_1_6=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\80_90_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
ObjV_7=ObjV;
parameter1_7=parameter1;
C_T_mean_1_7=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\100_grid_search_CSparameter1F1RIIN_0.1_0.1_2_EGFRICON_0.01_0.01_0.05_zhouqiwei100_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_8=ObjV;
parameter1_8=parameter1;
C_T_mean_1_8=C_T_mean_1;
ObjV_11=zeros(20,10);
ObjV_11(:,1:2)=ObjV_1(:,1:2);
ObjV_11(:,3)=ObjV_1(:,1);
ObjV_11(:,4)=ObjV_1(:,1);
ObjV_11(:,5)=ObjV_1(:,1);
ObjV_11(:,6)=ObjV_1(:,1);
ObjV_11(:,7)=ObjV_1(:,1);
ObjV_11(:,8:9)=ObjV_1(:,1:2);
ObjV_11(:,10)=ObjV_1(:,1);
parameter1_11=zeros(3,3);
parameter1_11(1,1)=0.1;parameter1_11(1,2)=0.1;parameter1_11(1,3)=2;
parameter1_11(2,1)=0.01;parameter1_11(2,2)=0.01;parameter1_11(2,3)=0.01;
parameter1_11(3,1)=10;parameter1_11(3,2)=10;parameter1_11(3,3)=100;
C_T_mean_1_11=zeros(301,20,10);
C_T_mean_1_11(:,:,1:2)=C_T_mean_1_1(:,:,1:2);
C_T_mean_1_11(:,:,3)=C_T_mean_1_2(:,:,1);
C_T_mean_1_11(:,:,4)=C_T_mean_1_3(:,:,1);
C_T_mean_1_11(:,:,5)=C_T_mean_1_4(:,:,1);
C_T_mean_1_11(:,:,6)=C_T_mean_1_5(:,:,1);
C_T_mean_1_11(:,:,7)=C_T_mean_1_6(:,:,1);
C_T_mean_1_11(:,:,8:9)=C_T_mean_1_7(:,:,1:2);
C_T_mean_1_11(:,:,10)=C_T_mean_1_8(:,:,1);
ObjV=ObjV_11;
parameter1=parameter1_11;
C_T_mean_1=C_T_mean_1_11;

% save('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\10_10_100_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');


load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\10_20_grid_search_CSF1RIIN_0.1_0.1_2_IGF1RICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_1=ObjV;
parameter1_1=parameter1;
C_T_mean_1_1=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\30_grid_search_CSF1RIIN_0.1_0.1_2_IGFRICON_0.01_0.01_0.05_zhouqiwei30_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_2=ObjV;
parameter1_2=parameter1;
C_T_mean_1_2=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\40_grid_search_CSF1RIIN_0.1_0.1_2_IGF1RICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_3=ObjV;
parameter1_3=parameter1;
C_T_mean_1_3=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\50_grid_search_CSF1RIIN_0.1_0.1_2_IGFRICON_0.01_0.01_0.05_zhouqiwei50_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_4=ObjV;
parameter1_4=parameter1;
C_T_mean_1_4=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\60_grid_search_CSF1RIIN_0.1_0.1_2_IGF1RICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
ObjV_5=ObjV;
parameter1_5=parameter1;
C_T_mean_1_5=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\70_grid_search_CSparameter1F1RIIN_0.1_0.1_2_IGFRICON_0.01_0.01_0.05_zhouqiwei70_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_6=ObjV;
parameter1_6=parameter1;
C_T_mean_1_6=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\80_90_grid_search_CSF1RIIN_0.1_0.1_2_IGF1RICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
ObjV_7=ObjV;
parameter1_7=parameter1;
C_T_mean_1_7=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\100_grid_search_CSparameter1F1RIIN_0.1_0.1_2_IGFRICON_0.01_0.01_0.05_zhouqiwei100_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_8=ObjV;
parameter1_8=parameter1;
C_T_mean_1_8=C_T_mean_1;

ObjV_11=zeros(20,10);
ObjV_11(:,1:2)=ObjV_1(:,1:2);
ObjV_11(:,3)=ObjV_1(:,1);
ObjV_11(:,4)=ObjV_1(:,1);
ObjV_11(:,5)=ObjV_1(:,1);
ObjV_11(:,6)=ObjV_1(:,1);
ObjV_11(:,7)=ObjV_1(:,1);
ObjV_11(:,8:9)=ObjV_1(:,1:2);
ObjV_11(:,10)=ObjV_1(:,1);
parameter1_11=zeros(3,3);
parameter1_11(1,1)=0.1;parameter1_11(1,2)=0.1;parameter1_11(1,3)=2;
parameter1_11(2,1)=0.01;parameter1_11(2,2)=0.01;parameter1_11(2,3)=0.01;
parameter1_11(3,1)=10;parameter1_11(3,2)=10;parameter1_11(3,3)=100;
C_T_mean_1_11=zeros(301,20,10);
C_T_mean_1_11(:,:,1:2)=C_T_mean_1_1(:,:,1:2);
C_T_mean_1_11(:,:,3)=C_T_mean_1_2(:,:,1);
C_T_mean_1_11(:,:,4)=C_T_mean_1_3(:,:,1);
C_T_mean_1_11(:,:,5)=C_T_mean_1_4(:,:,1);
C_T_mean_1_11(:,:,6)=C_T_mean_1_5(:,:,1);
C_T_mean_1_11(:,:,7)=C_T_mean_1_6(:,:,1);
C_T_mean_1_11(:,:,8:9)=C_T_mean_1_7(:,:,1:2);
C_T_mean_1_11(:,:,10)=C_T_mean_1_8(:,:,1);

ObjV=ObjV_11;
parameter1=parameter1_11;
C_T_mean_1=C_T_mean_1_11;
% save('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\10_10_100_grid_search_CSF1RIIN_0.1_0.1_2_IGFRICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');


% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_1.8_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_9=ObjV;
% parameter1_9=parameter1;
% C_T_mean_1_9=C_T_mean_1;
% load('D:\Matlab代码存放\2018胶质瘤\总程序\服务器数据2\grid_search_CSF1RIIN_2_IGF1RICON_0.2-0.2-2_Jian_zhouqi3.mat','ObjV','parameter1','C_T_mean_1');
% ObjV_10=ObjV;
% parameter1_10=parameter1;
% C_T_mean_1_10=C_T_mean_1;

load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.01_0.01_0.08_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_1=ObjV;
parameter1_1=parameter1;
C_T_mean_1_1=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.09_0.01_0.15_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
ObjV_2=ObjV;
parameter1_2=parameter1;
C_T_mean_1_2=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.16_0.01_0.2_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_3=ObjV;
parameter1_3=parameter1;
C_T_mean_1_3=C_T_mean_1;

C_T_mean_1_11=zeros(301,20,10);
C_T_mean_1_11(:,1:8,:)=C_T_mean_1_1(:,:,:);
C_T_mean_1_11(:,9:15,:)=C_T_mean_1_2(:,:,:);
C_T_mean_1_11(:,16:20,:)=C_T_mean_1_3(:,:,:);


C_T_mean_1=C_T_mean_1_11;
% save('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.01_0.01_0.2_EGFRICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');

load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.01_0.01_0.08_IGF1RICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_1=ObjV;
parameter1_1=parameter1;
C_T_mean_1_1=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.09_0.01_0.15_IGF1RICON_0.01_0.01_0.02_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_2=ObjV;
parameter1_2=parameter1;
C_T_mean_1_2=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.16_0.01_0.2_IGF1RICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
ObjV_3=ObjV;
parameter1_3=parameter1;
C_T_mean_1_3=C_T_mean_1;

C_T_mean_1_11=zeros(301,20,10);
C_T_mean_1_11(:,1:8,:)=C_T_mean_1_1(:,:,:);
C_T_mean_1_11(:,9:15,:)=C_T_mean_1_2(:,:,1,:);
C_T_mean_1_11(:,16:20,:)=C_T_mean_1_3(:,:,:);


C_T_mean_1=C_T_mean_1_11;
% save('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.01_0.01_0.2_IGF1RICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');

% ObjV=zeros(10,10);
% parameter1=zeros(2,3);
% C_T_mean_1=zeros(301,10,10);

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


%画图1
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\grid_search_CSF1RICON_0.1_0.1_2_EGFRICON_0.01.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_111=C_T_mean_1;
C_T_mean_11=C_T_mean_111(:,1:10);
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\10_10_100_grid_search_CSF1RIIN_0.1_0.1_2_EGFRICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_121=C_T_mean_1;
C_T_mean_12=C_T_mean_121(:,1:10,:);
C_T_mean_1=zeros(301,10,11);
C_T_mean_1(:,:,1)=C_T_mean_11;
C_T_mean_1(:,:,2:11)=C_T_mean_12;
C_T_mean_1_min_value=zeros(10,11);
C_T_mean_1_min_position=zeros(10,11);
C_T_mean_1_end_value=zeros(10,11);
for j=1:1:10
    for jj=1:1:11
        [aa1,bb1]=min(C_T_mean_1(:,j,jj));
        C_T_mean_1_min_value(j,jj)=aa1;
        C_T_mean_1_min_position(j,jj)=bb1;
        C_T_mean_1_end_value(j,jj)=C_T_mean_1(301,j,jj);
    end
end
C_T_mean_1_reappear=301*ones(10,11);
for j=1:1:10
    for jj=1:1:11
        cc1=C_T_mean_1_min_position(j,jj);
        cc2=find(C_T_mean_1(cc1:301,j,jj)>=0.6);
        if isempty(cc2)
            cccc=0
        else
            C_T_mean_1_reappear(j,jj)=C_T_mean_1_min_position(j,jj)+cc2(1)-1;
        end
    end
end
for j=1:1:10
    for jj=1:1:11      
        if  C_T_mean_1_reappear(j,jj)==1
            C_T_mean_1_reappear(j,jj)=301
        end
    end
end

[X1,Y1]=meshgrid(0:10:100,0.1:0.1:1);
[X2,Y2]=meshgrid(0:1:100,0.1:0.01:1);
ZI=interp2(X1,Y1,C_T_mean_1_reappear,X2,Y2,'cubic');


figure
 pcolor(1:101,1:91,ZI)
 title('CSF1RI & EGFRI Dose 0.01');
% colorbar('ylim',[0.7,0.78],'ytick',(0.7:0.01:0.78)); 
colorbar;
shading interp;
   set(gca,'ytick',1:10:91); set(gca,'yticklabel',0.1:0.1:1);
   set(gca,'xtick',1:20:101); set(gca,'xticklabel',0:20:100);
   xlabel('Period');,ylabel('CSF1RI');
   set(gca,'FontSize',14);

%画图2
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\grid_search_CSF1RICON_0.1_0.1_2_IGF1RICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_111=C_T_mean_1;
C_T_mean_11=C_T_mean_111(:,1:10);
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_20_30_100\10_10_100_grid_search_CSF1RIIN_0.1_0.1_2_IGFRICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_121=C_T_mean_1;
C_T_mean_12=C_T_mean_121(:,1:10,:);
C_T_mean_1=zeros(301,10,11);
C_T_mean_1(:,:,1)=C_T_mean_11;
C_T_mean_1(:,:,2:11)=C_T_mean_12;
C_T_mean_1_min_value=zeros(10,11);
C_T_mean_1_min_position=zeros(10,11);
C_T_mean_1_end_value=zeros(10,11);
for j=1:1:10
    for jj=1:1:11
        [aa1,bb1]=min(C_T_mean_1(:,j,jj));
        C_T_mean_1_min_value(j,jj)=aa1;
        C_T_mean_1_min_position(j,jj)=bb1;
        C_T_mean_1_end_value(j,jj)=C_T_mean_1(301,j,jj);
    end
end
C_T_mean_1_reappear=301*ones(10,11);
for j=1:1:10
    for jj=1:1:11
        cc1=C_T_mean_1_min_position(j,jj);
        cc2=find(C_T_mean_1(cc1:301,j,jj)>=0.6);
        if isempty(cc2)
            cccc=0
        else
            C_T_mean_1_reappear(j,jj)=C_T_mean_1_min_position(j,jj)+cc2(1)-1;
        end
    end
end
for j=1:1:10
    for jj=1:1:11      
        if  C_T_mean_1_reappear(j,jj)==1
            C_T_mean_1_reappear(j,jj)=301
        end
    end
end

[X1,Y1]=meshgrid(0:10:100,0.1:0.1:1);
[X2,Y2]=meshgrid(0:1:100,0.1:0.01:1);
ZI=interp2(X1,Y1,C_T_mean_1_reappear,X2,Y2,'cubic');

figure
 pcolor(1:101,1:91,ZI)
 title('CSF1RI & IGF1RI Dose 0.01');
% colorbar('ylim',[0.7,0.78],'ytick',(0.7:0.01:0.78)); 
colorbar;
shading interp;
   set(gca,'ytick',1:10:91); set(gca,'yticklabel',0.1:0.1:1);
   set(gca,'xtick',1:20:101); set(gca,'xticklabel',0:20:100);
   xlabel('Period');,ylabel('CSF1RI');
   set(gca,'FontSize',14);


%画图3
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\new_grid_search_CSF1RICON_0.01_0.01_0.2_EGFRICON_0.01_Jian.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_11=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.01_0.01_0.2_EGFRICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_12=C_T_mean_1;
C_T_mean_1=zeros(301,20,11);
C_T_mean_1(:,:,1)=C_T_mean_11(:,:,1);
C_T_mean_1(:,:,2:11)=C_T_mean_12;
C_T_mean_1_min_value=zeros(20,11);
C_T_mean_1_min_position=zeros(20,11);
C_T_mean_1_end_value=zeros(20,11);
for j=1:1:20
    for jj=1:1:11
        [aa1,bb1]=min(C_T_mean_1(:,j,jj));
        C_T_mean_1_min_value(j,jj)=aa1;
        C_T_mean_1_min_position(j,jj)=bb1;
        C_T_mean_1_end_value(j,jj)=C_T_mean_1(301,j,jj);
    end
end
C_T_mean_1_reappear=301*ones(20,11);
for j=1:1:20
    for jj=1:1:11
        cc1=C_T_mean_1_min_position(j,jj);
        cc2=find(C_T_mean_1(cc1:301,j,jj)>=0.6);
        if isempty(cc2)
            cccc=0
        else
            C_T_mean_1_reappear(j,jj)=C_T_mean_1_min_position(j,jj)+cc2(1)-1;
        end
    end
end
for j=1:1:20
    for jj=1:1:11      
        if  C_T_mean_1_reappear(j,jj)==1
            C_T_mean_1_reappear(j,jj)=301
        end
    end
end

[X1,Y1]=meshgrid(0:10:100,0.01:0.01:0.2);
[X2,Y2]=meshgrid(0:1:100,0.01:0.001:0.2);
ZI=interp2(X1,Y1,C_T_mean_1_reappear,X2,Y2,'cubic');


figure
 pcolor(1:101,1:191,ZI)
 title('CSF1RI & EGFRI Dose 0.01');
% colorbar('ylim',[0.7,0.78],'ytick',(0.7:0.01:0.78)); 
colorbar;
shading interp;
   set(gca,'ytick',1:10:191); set(gca,'yticklabel',0.01:0.01:0.2);
   set(gca,'xtick',1:20:101); set(gca,'xticklabel',0:20:100);
   xlabel('Period');,ylabel('CSF1RI');
   set(gca,'FontSize',14);

%画图4
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\grid_search_CSF1RICON_0.01_0.01_0.2_IGF1RICON_0.01_0.01_0.02.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_11=C_T_mean_1;
load('D:\Matlab代码存放\2018胶质瘤\总程序\11月28\10_10_100\10_10_100_grid_search_CSF1RIIN_0.01_0.01_0.2_IGF1RICON_0.01_jian.mat','ObjV','parameter1','C_T_mean_1');
C_T_mean_12=C_T_mean_1;
C_T_mean_1=zeros(301,20,11);
C_T_mean_1(:,:,1)=C_T_mean_11(:,:,1);
C_T_mean_1(:,:,2:11)=C_T_mean_12;
C_T_mean_1_min_value=zeros(20,11);
C_T_mean_1_min_position=zeros(20,11);
C_T_mean_1_end_value=zeros(20,11);
for j=1:1:20
    for jj=1:1:11
        [aa1,bb1]=min(C_T_mean_1(:,j,jj));
        C_T_mean_1_min_value(j,jj)=aa1;
        C_T_mean_1_min_position(j,jj)=bb1;
        C_T_mean_1_end_value(j,jj)=C_T_mean_1(301,j,jj);
    end
end
C_T_mean_1_reappear=301*ones(20,11);
for j=1:1:20
    for jj=1:1:11
        cc1=C_T_mean_1_min_position(j,jj);
        cc2=find(C_T_mean_1(cc1:301,j,jj)>=0.6);
        if isempty(cc2)
            cccc=0
        else
            C_T_mean_1_reappear(j,jj)=C_T_mean_1_min_position(j,jj)+cc2(1)-1;
        end
    end
end
for j=1:1:20
    for jj=1:1:11      
        if  C_T_mean_1_reappear(j,jj)==1
            C_T_mean_1_reappear(j,jj)=301
        end
    end
end

[X1,Y1]=meshgrid(0:10:100,0.01:0.01:0.2);
[X2,Y2]=meshgrid(0:1:100,0.01:0.001:0.2);
ZI=interp2(X1,Y1,C_T_mean_1_reappear,X2,Y2,'cubic');


figure
 pcolor(1:101,1:191,ZI)
 title('CSF1RI & IGF1RI Dose 0.01');
% colorbar('ylim',[0.7,0.78],'ytick',(0.7:0.01:0.78)); 
colorbar;
shading interp;
   set(gca,'ytick',1:10:191); set(gca,'yticklabel',0.01:0.01:0.2);
   set(gca,'xtick',1:20:101); set(gca,'xticklabel',0:20:100);
   xlabel('Period');,ylabel('CSF1RI');
   set(gca,'FontSize',14);



C_T_mean_1_reappear_2=zeros(10,10);
for j=1:1:10
    for jj=1:1:10
        cc1=C_T_mean_1_min_position(j,jj);
        cc2=find(C_T_mean_1(cc1:301,j,jj)>=0.6);
        if isempty(cc2)
            cccc=0
        else
            C_T_mean_1_reappear_2(j,jj)=cc2(1);
        end
    end
end