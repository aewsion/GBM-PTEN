%%% LHS_PRCC
clear all
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
global eta_d D_d
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
global dmax1 dmax2 dmax3
global select_CSF1R_I select_EGFR_I select_IGF1R_I a_A
global b_EGF a_M2M1

Parameters_nondimensional
PTEN=0;
AKT_I=0;
GSK_I=0;

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters
Parameter=[alpha_AKT K5 B_Gal9 K2 d_TM1 B_M1M2];

%%%%  LHS

p = 1000;   % Number of points
num_para = length(Parameter);   % Number of dimensions
lb = Parameter*0.9; % lower bounds for parameters
ub = Parameter*1.1; % upper bounds for parameters
X = lhsdesign(p,num_para,'criterion','correlation');
D = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb)));%plus应该是相加的意思，times应该是相乘的意思

r_RG=zeros(1,p);
r_K=zeros(1,p);

C=zeros(p,TimeLength+1);
%%%%%%%  Model Simulation
for i=1:1:p
    i
    % Parameter=[alpha_AKT K5 B_Gal9 K2 A_Gal9 d_TM1 B_M1M2];
    alpha_AKT=D(i,1);
    K5=D(i,2);
    B_Gal9=D(i,3);
    K2=D(i,4);
    d_TM1=D(i,5);
    B_M1M2=D(i,6);


    options=odeset('reltol',1e-3);%精度函数不能少
    m = 2;%用matlab解需要设置的参数
    sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);
    u1 = sol(:,:,1);%C_T
    C_T_mean=sum(u1(:,:),2)/(L2+1); %0-200day的肿瘤平均密度
    C(i,:)=C_T_mean;%1000*200
    timelength=size(C,2);



    if C(i,timelength)~=min(C(i,:))
        r_RG(i)=(C(i,timelength)-min(C(i,:)))/(timelength-find(C(i,:)==min(C(i,:)),1));
    else
        r_RG(i)=0;
    end
    if C(i,1)~=min(C(i,:))
        r_K(i)=(C(i,1)-min(C(i,:)))/(find(C(i,:)==min(C(i,:)),1));
    else
        r_K(i)=0;%平均每天肿瘤减少多少？
    end
    C_negative=(C<0);

end

[asizeC,bsizeC]=size(C);
TimeLength=bsizeC-1;
threshold=0.65;
[p,p12]=size(C);
death_day=zeros(p,1);
%death_day表示死亡的天数，death_person表示每天死亡的人数，death_person_20表示每20天死亡的人数
for person=1:1:p
    if isempty(find(C(person,:)>=threshold,1))%为空进入
        max_death_day=max(C(person,:));
        i=find(C(person,:)==max_death_day);
        if i+50<=p12 && max_death_day>=0.6 && (max_death_day-sum(C(person,i+1:i+50))/50)/max_death_day<0.005
            death_day(person,1)=i;
        end
    else
        death_day(person,1)=find(C(person,:)>=threshold,1); %只要第一个，大于阈值的第一天
    end
end
death_person=zeros(TimeLength+1,1);
for i=1:1:TimeLength+1
    dp=find(death_day==i);
    if isempty(dp)%空就是这天无人死亡
        death_person(i,1)=length(dp);
    end
end
Survival_Rate=100*ones(TimeLength+1,1);
%death_day表示死亡的天数，death_person表示每天死亡的人数，death_person_20表示每20天死亡的人数
for i=2:1:TimeLength+1
    Survival_Rate(i,1)=(p-death_person(i,1))*100/p;
end
death_person_20=zeros(TimeLength/20+1,1);
for i=1:1:TimeLength/20
    death_person_20(i+1,1)=sum(death_person((i-1)*20+1:i*20,1));
end

save('KO_1000.mat','C','r_RG','r_K','D','death_day','death_person','death_person_20');


load('KO_1000.mat','C','r_RG','r_K','D','death_day','death_person','death_person_20');
% load('wild_1000.mat','C','r_RG','r_K','D','death_day','death_person','death_person_20');

[p,num_para]=size(D);%1000*参数个数
if num_para==6
    %%% Ranking
    alpha_AKT=D(:,1);
    K5=D(:,2);
    B_Gal9=D(:,3);
    K2=D(:,4);
    d_TM1=D(:,5);
    B_M1M2=D(:,6);
    parameter=[alpha_AKT,K5,B_Gal9,K2,d_TM1,B_M1M2];
    Para_label={'{\alpha}_{AKT}','K5','B_{Gal9}','K2','d_{TM}','B_{M1M2}'};
elseif num_para==11
    alpha_AKT=D(:,1);
    K5=D(:,2);
    B_Gal9=D(:,3);
    K2=D(:,4);
    d_TM1=D(:,5);
    B_M1M2=D(:,6);
    parameter=[alpha_AKT,K5,B_Gal9,K2,A_Gal9,d_TM1,B_M1M2];
    Para_label={'{\alpha}_{AKT}','K5','B_{Gal9}','K2','d_{TM}','B_{M1M2}'};
else
    error('Cannot calculate with given values')
end

r_RG_O=r_RG';%转置
r_K_O=r_K';

[P_R,ind_P]=sort(parameter);%sort对每列排序,参数一排序后矩阵，参数二排序后所对应的原坐标矩阵
[r_RG_R,ind_RG]=sort(r_RG_O);
[r_K_R,ind_K]=sort(r_K_O);


[r_RG_R,ind_RG]=sort(death_day);%死亡天数


CC=zeros(num_para,1);
Pvalue=zeros(num_para,1);
%%%%  Pearson CC for raw data
for j=1:num_para%num_para参数个数
    a=1:num_para;
    a(a==j)=[];
    %         [rho,pval] = corr(r_RG_O,parameter(:,j));
    [rho,pval] = corr(death_day,parameter(:,j));
    CC(j)=rho;
    Pvalue(j)=pval;
end

b1=[110/255,165/255,211/255];b2=[241/255,153/255,157/255];
figure,
bar(1:num_para,CC)  %相关系数
set(gca,'xtick',1:1:num_para)
set(gca,'xticklabel',Para_label)
rotateticklabel(gca,315);
% xlabel('Parameters','FontWeight','Bold','FontSize',11);
ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
title('Correlation Coefficient','FontWeight','Bold','FontSize',11);
set(gca,'FontSize',14);
figure_FontSize=11;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',14),'FontSize',figure_FontSize);

figure,
bar(1:num_para,Pvalue)  %
set(gca,'xtick',1:1:num_para)
set(gca,'xticklabel',Para_label)
rotateticklabel(gca,315);
% xlabel('Parameters','FontWeight','Bold','FontSize',11);
% ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
title('Correlation Coefficient Pvalue','FontWeight','Bold','FontSize',11);
set(gca,'FontSize',14);
figure_FontSize=11;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',14),'FontSize',figure_FontSize);



% %%%%  Partial CC for raw data
% %[R,P] = partialcorr(X,Y,Z); %在控制变量Z的影响下，计算变量X、Y的偏相关系数。
% for j=1:num_para
%     a=1:num_para;
%     a(a==j)=[];
%     [rho_RG,pval_RG] = partialcorr(r_RG_O,parameter(:,j),parameter(:,a));
%     [rho_K,pval_K] = partialcorr(r_K_O(~isnan(r_K_O)),parameter(~isnan(r_K_O),j),parameter(~isnan(r_K_O),a));
%     PCC_RG(j)=rho_RG;
%     PCC_K(j)=rho_K;
%     Pvalue_RG(j)=pval_RG;
%     Pvalue_K(j)=pval_K;
% end
% figure,
% hh1=bar(1:num_para,PCC_RG);  %
% set(hh1,'facecolor',b1);
% set(gca,'xtick',1:1:num_para)
% set(gca,'xticklabel',Para_label)
% rotateticklabel(gca,315);
% % xlabel('Parameters','FontWeight','Bold','FontSize',11);
% ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
% title('Regrowth Rate','FontWeight','Bold','FontSize',11);
% set(gca,'FontSize',11);
% figure_FontSize=11;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',14),'FontSize',figure_FontSize);
% 
% figure,
% bar(1:num_para,Pvalue_RG);  %
% set(gca,'xtick',1:1:num_para)
% set(gca,'xticklabel',Para_label)
% rotateticklabel(gca,315);
% % xlabel('Parameters','FontWeight','Bold','FontSize',11);
% ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
% title('Pvalue RG','FontWeight','Bold','FontSize',11);
% set(gca,'FontSize',14);
% figure_FontSize=11;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',14),'FontSize',figure_FontSize);
% 
% figure,
% hh2=bar(1:num_para,PCC_K);  %
% set(hh2,'facecolor',b2);
% set(gca,'xtick',1:1:num_para)
% set(gca,'xticklabel',Para_label)
% rotateticklabel(gca,315);
% % xlabel('Parameters','FontWeight','Bold','FontSize',11);
% ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
% title('Killing Rate','FontWeight','Bold','FontSize',11);
% set(gca,'FontSize',11);
% figure_FontSize=11;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',14),'FontSize',figure_FontSize);
% 
% figure,
% bar(1:num_para,Pvalue_K)  %
% set(gca,'xtick',1:1:num_para)
% set(gca,'xticklabel',Para_label)
% rotateticklabel(gca,315);
% % xlabel('Parameters','FontWeight','Bold','FontSize',11);
% ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
% title('Pvalue K','FontWeight','Bold','FontSize',11);
% set(gca,'FontSize',14);
% figure_FontSize=11;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',14),'FontSize',figure_FontSize);
% 
% %%%%  Partial CC for ranks
% for j=1:num_para
%     a=1:num_para;
%     a(a==j)=[];
%     [rho,pval] = partialcorr(ind_RG,ind_P(:,j),ind_P(:,a));
%     PRCC(j)=rho;
%     Pvalue(j)=pval;
% end
% figure,
% bar(1:num_para,PRCC)  %
% set(gca,'xtick',1:1:num_para)
% set(gca,'xticklabel',Para_label)
% rotateticklabel(gca,315);
% % xlabel('Parameters','FontWeight','Bold','FontSize',11);
% ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
% title('Partial Correlation Coefficients For Rank','FontWeight','Bold','FontSize',11);
% set(gca,'FontSize',14);
% figure_FontSize=11;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',14),'FontSize',figure_FontSize);
% 
% figure,
% bar(1:num_para,Pvalue)  %
% set(gca,'xtick',1:1:num_para)
% set(gca,'xticklabel',Para_label)
% rotateticklabel(gca,315);
% % xlabel('Parameters','FontWeight','Bold','FontSize',11);
% ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',11);
% title('Partial Correlation Coefficients Pvalue For Rank','FontWeight','Bold','FontSize',11);
% set(gca,'FontSize',14);
% figure_FontSize=11;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',14),'FontSize',figure_FontSize);


toc
