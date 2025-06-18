% === This test program is design for estimating ODE system parameters in Signaling Pathway====
% === based on RPPA data.              ====
%
% Written by Xiaoqiang Sun   Apr 05, 2013
% Wake Forerst Medicel School, Zhou's Lab 
% Winston Salem, NC, U.S.A 27157 
tic

clc
clear all
close all
figure,
global time
time=0;
global Ptest1 

group=10;%想要用遗传算法来估计10组参数
NVAR =22;%model 2            % Number of variables, according to the original setting of Para0  每组参数一共有22个参数
p=ones(group,NVAR+1);     %p来保存10组参数的22个参数加上每一组参数的评估值（计算得到），一共10×23大小
for q=1:1:group
NIND = 20;           % Number of individuals per subpopulations    相当于每次20条染色体
MAXGEN = 400;        % maximum Number of generations     每找一组参数迭代的次数 400次
GGAP = 0.90;         % Generation gap, how many new individuals are created
PRECI = 8;           % Precision of binary representation  精确到8位二进制
Para0= ones(1,NVAR);  % initial values in parameter space
Ptest1=Para0;
shu_min=0.000000001;           %Minimum search range    首先设置每个参数的最小值
shu_max=50;          %Maximum search range        再设置有的参数的最大值
% Build field descriptor
% FieldD = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,20) 1;shu_max*ones(1,10) 1*ones(1,8) 10*ones(1,2) 11]
%               rep([1; 0; 1 ;1], [1, NVAR])];%model 1
FieldD = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,21) 9;shu_max*ones(1,10) 1*ones(1,9) 10*ones(1,2) 10.2]
              rep([1; 0; 1 ;1], [1, NVAR])];%model 2    解码函数因为要二进制和十进制来转换运算
% FieldD = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,19) 1;shu_max*ones(1,10) 1*ones(1,7) 10*ones(1,2) 11]
%               rep([1; 0; 1 ;1], [1, NVAR])];%model 3
% FieldD = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,18) 1;shu_max*ones(1,10) 1*ones(1,6) 10*ones(1,2) 11]
%               rep([1; 0; 1 ;1], [1, NVAR])];%model 4

% Initialise population
Chrom = crtbp(NIND, NVAR*PRECI);   %8位二进制
% Reset counters
Best = NaN*ones(MAXGEN,1);	% best in current population
gen = 1;			        % generational counter
rvChrom=bs2rv(Chrom,FieldD);  %翻译染色体
ObjV = Objfun1(rvChrom);  %  %  %%% Objective function % Evaluate initial population
for i=1:1:MAXGEN          %遗传算法迭代400次寻找最优的染色体
    gen=i  %为了在命令框看当前迭代次数
     % Assign fitness-value to entire population
       FitnV = ranking(ObjV);
    % Select individuals for breeding   选择
       SelCh = select('sus', Chrom, FitnV, GGAP);
    % Recombine selected individuals (crossover)   交叉
       SelCh = recombin('xovsp',SelCh,0.7);
    % Perform mutation on offspring  变异
       SelCh = mut(SelCh);
    % Evaluate offspring, call objective function   评估每条染色体的计算得到的拟合评分
       ObjVSel = Objfun1(bs2rv(SelCh,FieldD));
    % Reinsert offspring into current population  得到新的群体也就是新的20条染色体
       [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);


    % Update display and record current best individual  找出最优的那条染色体
       [Best(i),ind] = min(ObjV);
      %%画图
       plot(log10(Best),'ro'); xlabel('generation'); ylabel('log10(f(x))');
       text(0.5,0.95,['Best = ', num2str(Best(gen))],'Units','normalized');
       drawnow;     
        % Increment generational counter
end
chrom_value=bs2rv(Chrom,FieldD)   % Show the last set of subpopulation  翻译染色体
% for j=1:1:500
%下面是重新设置寻找参数的范围，最大值和最小值，这部分也可以不用，直接在前一步增加迭代次数也能搜索出来，更新搜索范围为了更好搜索
[value ind]=min(ObjV);
Ptest2=chrom_value(ind,:)         
for j=1:1:10
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,100);
end
for j=11:1:19
    FieldD(2,j)=max(Ptest2(j)-0.1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+0.1,1);     
end
for j=20:1:21
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,10);       
end
for j=22:1:22
    FieldD(2,j)=max(Ptest2(j)-2,1);
    FieldD(3,j)=min(Ptest2(j)+2,10.2);
end
b_min=FieldD(2,:);
b_max=FieldD(3,:);
%min_old+(max_old-min_old)*rand=min_new+(max_new-min_new)*rand_new
chrom_new=(Ptest2-b_min)./(b_max-b_min);
chrom_new2=chrom_new*2^8;
chrom_new2=floor(chrom_new2);
for i=0:1:7
    b=floor(chrom_new2/2^(7-i));
    for j=0:1:NVAR-1   
    Chrom(ind,((8-(7-i))+8*j))=b(j+1);
    end
    chrom_new2=rem(chrom_new2,2^(7-i));
end
rvChrom=bs2rv(Chrom,FieldD);
ObjV = Objfun1(rvChrom);  %  %  %%% Objective function % Evaluate initial population
%下面是更新搜索范围后开始小迭代次数迭代，然后再更新搜索范围，一共100次更新最后搜索出最优的染色体，下面也是遗传算法
for i=1:1:30
    %i为搜索的次数
    gen=1;
    MAXGEN=30;
    Best = NaN*ones(MAXGEN,1);	% best in current population
    rvChrom=bs2rv(Chrom,FieldD);
    ObjV = Objfun1(rvChrom);  %  %  %%% Objective function % Evaluate initial population
while gen < MAXGEN+1
    gen%方便看运行到哪里设置的
    % Assign fitness-value to entire population
       FitnV = ranking(ObjV);
    % Select individuals for breeding
       SelCh = select('sus', Chrom, FitnV, GGAP);
    % Recombine selected individuals (crossover)
       SelCh = recombin('xovsp',SelCh,0.7);
    % Perform mutation on offspring
       SelCh = mut(SelCh);
    % Evaluate offspring, call objective function
       ObjVSel = Objfun1(bs2rv(SelCh,FieldD));

    % Reinsert offspring into current population
    %将后代重新插入当前种群
       [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);


    % Update display and record current best individual
       [Best(gen),ind] = min(ObjV);
      
       
       plot(log10(Best),'ro'); xlabel('generation'); ylabel('log10(f(x))');
       text(0.5,0.95,['Best = ', num2str(Best(gen))],'Units','normalized');
       drawnow;%更新图窗并处理任何挂起的回调。如果您修改图形对象并且需要在屏幕上立即查看这次更新，请使用该命令。       
        % Increment generational counter
       gen = gen+1;
end 
chrom_value=bs2rv(Chrom,FieldD)   % Show the last set of subpopulation
% for j=1:1:500
[value,ind]=min(ObjV);
Ptest2=chrom_value(ind,:)         %把最优个体找出来
Ptest1=Ptest2
for j=1:1:10
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,100);
end
for j=11:1:19
    FieldD(2,j)=max(Ptest2(j)-0.1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+0.1,1);     
end
for j=20:1:21
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,10);       
end
for j=22:1:22
    FieldD(2,j)=max(Ptest2(j)-2,1);
    FieldD(3,j)=min(Ptest2(j)+2,10.2);
end
b_min=FieldD(2,:);
b_max=FieldD(3,:);
%min_old+(max_old-min_old)*rand=min_new+(max_new-min_new)*rand_new
chrom_new=(Ptest2-b_min)./(b_max-b_min);
chrom_new2=chrom_new*2^8;
chrom_new2=floor(chrom_new2);
for k=0:1:7
    b=floor(chrom_new2/2^(7-k));
    for j=0:1:NVAR-1   %NVAR为参数个数
    Chrom(ind,((8-(7-k))+8*j))=b(j+1);
    end
    chrom_new2=rem(chrom_new2,2^(7-k));
end
%    b_min=FieldD(2,:);
%    b_max=FieldD(3,:);
%    b_min2=b_min(ones(NIND,1),:);
%    b_max2=b_max(ones(NIND,1),:);
%    %min_old+(max_old-min_old)*rand=min_new+(max_new)*rand_new
%    chrom_value_new=(chrom_value-b_min2)./(b_max2-b_min2);
  
end
%  Ptest1=p(1,1:23)
%  Conditions_ODES_2;
%  tick={'Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906'};
%  figure,b=bar([Y2(:,1) pY1068_EGFR]);title('pY1068-EGFR');
%  set(b(1,1),'FaceColor',[242/255,212/255,201/255]);
%  set(b(1,2),'FaceColor',[177/255,86/255,70/255]);
%  set(gca,'XTickLabelRotation',45);
%  set(gca,'Xticklabel',tick);
% %  xlabel('Condition','FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  set(gca,'FontWeight','Bold','FontSize',12);
%  legend('Model','Experiments');
%  figure, b=bar([Y2(:,2) p_ERK]);title('p-ERK');
%  set(b(1,1),'FaceColor',[242/255,203/255,5/255]);
%  set(b(1,2),'FaceColor',[242/255,159/255,5/255]);
%  set(gca,'XTickLabelRotation',45);
%  set(gca,'Xticklabel',tick);
% %  xlabel('Condition','FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  set(gca,'FontWeight','Bold','FontSize',12);
%  legend('Model','Experiments');
%  figure, b=bar([Y2(:,3) p_IGF1Rbeta]);title('p-IGF1Rbeta');
%  set(b(1,1),'FaceColor',[242/255,212/255,201/255]);
%  set(b(1,2),'FaceColor',[177/255,86/255,70/255]);
%  set(gca,'XTickLabelRotation',45);
%  set(gca,'Xticklabel',tick);
% %  xlabel('Condition','FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  set(gca,'FontWeight','Bold','FontSize',12);
%  legend('Model','Experiments');
%  figure, b=bar([Y2(:,4) pS473_AKT]);title('pS473-AKT'); 
%  set(b(1,1),'FaceColor',[242/255,203/255,5/255]);
%  set(b(1,2),'FaceColor',[242/255,159/255,5/255]);
%  set(gca,'XTickLabelRotation',45);
%  set(gca,'Xticklabel',tick);
% %  xlabel('Condition','FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  set(gca,'FontWeight','Bold','FontSize',12);
%  legend('Model','Experiments');
 Conditions_ODES_2;
 tick={'Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906'};
 b1=[147/255,224/255,255/255];b2=[153/255,77/255,82/255];
 figure,plot(1:c2,Y2(:,1),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,pY1068_EGFR,'-o','color',b2,'MarkerFaceColor',b2);
  %set(gca,'ytick',[0.1:2]);
%  xlabel('Condition','FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
 set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure,plot(1:c2,Y2(:,2),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,p_ERK,'-o','color',b2,'MarkerFaceColor',b2);
 set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
 set(gca,'YTick',0.4:0.1:1); 
 set(gca,'Yticklabel',0.4:0.1:1,'FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('p-ERK','FontWeight','Bold','FontSize',18);
 set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure,plot(1:c2,Y2(:,3),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,p_IGF1Rbeta,'-o','color',b2,'MarkerFaceColor',b2);
 set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
%  xlabel('Condition','FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
 set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure,plot(1:c2,Y2(:,4),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,pS473_AKT,'-o','color',b2,'MarkerFaceColor',b2); 
  set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
%  xlabel('Condition','FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('pS473-AKT','FontWeight','Bold','FontSize',18);
  set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure
 p(q,1:NVAR)=Ptest2;
 p(q,NVAR+1)=value;
end
save('canshu1-24hours-ode15s_wode.mat','p');
 
% value= Objfun1(Ptest1)

% save Ptest2.mat  % save the estimated parameters   %Ptest2才是最终结果
% load 'Ptest2.mat'

% save('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终\canshu1.mat','p') 
% save('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终2\canshu119.mat','p');
% save('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终2\canshu114-24hours-ode15s.mat','p');
% save('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终\canshu3.mat','p')
% save('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终\canshu4.mat','p')
% load('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终2\canshu114-24hours-ode15s.mat','p');
% load('canshu1-24hours-ode15s.mat','p');
% load('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终\canshu3.mat','p');
% load('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终\canshu3.mat','p');
% p2=p(1:4,:)
% p(7,:)=p2;
% p3=zeros(12,23)
% p3(1:10,:)=p;
% p3(9:12,:)=p2
% p=p3

% % load('canshu1-24hours-ode15s.mat','p');
% %画图
% % % Ptest1=p(1,1:22);
% for i=1:1:3
%     Ptest1=p(i,:)
%     Conditions_ODES_2;
%  tick={'Control 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906'};
%  b1=[0/255,90/255,171/255];b2=[255/255,66/255,93/255];
%  figure,plot(1:c2,Y2(:,1),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,pY1068_EGFR,'-o','color',b2,'MarkerFaceColor',b2);
%   %set(gca,'ytick',[0.1:2]);
% %  xlabel('Condition','FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  set(gca,'XTick',1:1:11);
%  set(gca,'XTickLabelRotation',45);
%  legend('Model','Experiments');
%  set(gca,'FontWeight','Bold','FontSize',18);
%  title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
%  set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
%  figure,plot(1:c2,Y2(:,2),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,p_ERK,'-o','color',b2,'MarkerFaceColor',b2);
%  set(gca,'XTick',1:1:11);
%  set(gca,'XTickLabelRotation',45);
%  set(gca,'YTick',0.4:0.2:1); 
%  set(gca,'Yticklabel',0.4:0.2:1,'FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  legend('Model','Experiments');
%  set(gca,'FontWeight','Bold','FontSize',18);
%  title('p-ERK','FontWeight','Bold','FontSize',18);
%  set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
%  figure,plot(1:c2,Y2(:,3),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,p_IGF1Rbeta,'-o','color',b2,'MarkerFaceColor',b2);
%  set(gca,'XTick',1:1:11);
%  set(gca,'XTickLabelRotation',45);
% %  xlabel('Condition','FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  legend('Model','Experiments');
%  set(gca,'FontWeight','Bold','FontSize',18);
%  title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
%  set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
%  figure,plot(1:c2,Y2(:,4),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,pS473_AKT,'-o','color',b2,'MarkerFaceColor',b2); 
%   set(gca,'XTick',1:1:11);
%  set(gca,'XTickLabelRotation',45);
% %  xlabel('Condition','FontWeight','Bold','FontSize',18);
%  ylabel('Relative value','FontWeight','Bold','FontSize',18);
%  legend('Model','Experiments');
%  set(gca,'FontWeight','Bold','FontSize',18);
%  title('pS473-AKT','FontWeight','Bold','FontSize',18);
%   set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
% end
% figure
% plot(0:1:21,p(1,1:22),'-+','LineWidth',1.5);hold on;
% plot(0:1:21,p(2,1:22),'--o','LineWidth',1.5);hold on;
% plot(0:1:21,p(3,1:22),':*','LineWidth',1.5);hold on;
% plot(0:1:21,p(4,1:22),'-.p','LineWidth',1.5);hold on;
% plot(0:1:21,p(5,1:22),'-x','LineWidth',1.5);hold on;
% plot(0:1:21,p(6,1:22),'--d','LineWidth',1.5);hold on;
% plot(0:1:21,p(7,1:22),':^','LineWidth',1.5);hold on;
% plot(0:1:21,p(8,1:22),'-.>','LineWidth',1.5);hold on;
% plot(0:1:21,p(9,1:22),'-h','LineWidth',1.5);hold on;
% plot(0:1:21,p(10,1:22),'--s','LineWidth',1.5);hold on;
%  title('Parameter 1-22 and the Objective function value','FontWeight','Bold','FontSize',18);
%  tick4={'V1','a','V21','V22','V3','V4','d1','d2','d3','d4','K11','K12','K21','K22','K23','K31','K32','K41','K42','K13','K33','n'};
%  set(gca,'XTickLabelRotation',45);
%  set(gca,'FontWeight','Bold','FontSize',18);
%  set(gca,'XTick',0:1:21);
%  set(gca,'Xticklabel',tick4,'FontWeight','Bold','FontSize',10);
%  set(gca,'XLim',[0 21]);
%  legend
%  
 
%  %敏感图
% %  load('canshu1-24hours-ode15s.mat','p');
% %  global Ptest1
% %  Ptest2=p(1,1:22);
% %  Ptest1=p(1,1:22);
% %  Q=zeros(11,25,4);
% %  Conditions_ODES_2;
% %  Q1=Q;
% %  q=zeros(22,4);
% %  for i=1:1:22
% %      Ptest1=Ptest2;
% %      Ptest1(1,i)=Ptest2(1,i)*1.01;
% %      Conditions_ODES_2;
% %       for j1=1:1:4
% %           for j2=1:1:11
% %              for j3=1:1:24
% %                  q(i,j1)=q(i,j1)+abs(Q1(j2,j3,j1)-Q(j2,j3,j1));
% %              end
% %           end
% %       end
% %  end
% %  q2=q/25/11/0.01;
% %  q3=q2*100;
% %  figure
% %  h=bar3(q3);
% %  colorbar; shading interp;
% %  tick1={'V1','a','V21','V22','V3','V4','d1','d2','d3','d4','K11','K12','K21','K22','K23','K31','K32','K41','K42','K13','K33'};
% %  tick2={'EGFR','ERK','IGF1R','AKT'};
% % %  set(gca,'YTickLabelRotation',15);
% % %  set(gca,'XTickLabelRotation',-10);
% %  set(gca,'YTick',1:1:21);
% %  set(gca,'Yticklabel',tick1);
% %  set(gca,'XTick',1:1:4);
% %  set(gca,'Xticklabel',tick2);
% %  set(gca,'FontWeight','Bold','FontSize',10);
% % %  ylabel('Parameters','FontWeight','Bold','FontSize',15,'rotation',-38);
% % %  xlabel('Variables','FontWeight','Bold','FontSize',12,'rotation',28);
% %  zlabel('Sensitivity(%)','FontWeight','Bold','FontSize',16);
% %  set(gca,'YLim',[0.5 21.5]);
% %  set(gcf,'color','w');
% %  for n=1:1:numel(h)
% %      cdata=get(h(n),'zdata');
% %      set(h(n),'cdata',cdata,'facecolor','interp')
% %  end
% % %  set(gca,'Position',get(gca,'OuterPosition')-get(gca,'TightInset')*[-1 0 1 0;0 -1 0 1;0 0 1 0;0 0 0 1]);
% % %  delete(gca);clf;
% %  Ptest2=p(2,1:22);
% %  Ptest1=Ptest2;
% %  Conditions_ODES_2;
% %  
% %   figure
% %   plot(0:1:48,Q(1,:,1),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(2,:,1),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(3,:,1),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(4,:,1),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(5,:,1),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(6,:,1),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(7,:,1),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(8,:,1),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(9,:,1),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(10,:,1),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(11,:,1),':','LineWidth',5);hold on;
% %   title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
% %   set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %   xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %   ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% %   legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% %   set(gca,'FontWeight','Bold','FontSize',18);
% %   set(gca,'XLim',[0 48]);
% %   
% %   figure
% %   plot(0:1:48,Q(1,:,2),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(2,:,2),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(3,:,2),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(4,:,2),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(5,:,2),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(6,:,2),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(7,:,2),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(8,:,2),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(9,:,2),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(10,:,2),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(11,:,2),':','LineWidth',5);hold on;
% %   title('p-ERK','FontWeight','Bold','FontSize',18);
% %   set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %   xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %   ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% %   legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% %   set(gca,'FontWeight','Bold','FontSize',18);
% %   set(gca,'XLim',[0 48]);
% %   
% %   figure
% %   plot(0:1:48,Q(1,:,3),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(2,:,3),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(3,:,3),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(4,:,3),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(5,:,3),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(6,:,3),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(7,:,3),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(8,:,3),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(9,:,3),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(10,:,3),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(11,:,3),':','LineWidth',5);hold on;
% %   title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
% %   set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %   xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %   ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% %   legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% %   set(gca,'FontWeight','Bold','FontSize',18);
% %   set(gca,'XLim',[0 48]);
% %   
% %   figure
% %   plot(0:1:48,Q(1,:,4),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(2,:,4),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(3,:,4),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(4,:,4),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(5,:,4),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(6,:,4),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(7,:,4),':','LineWidth',5);hold on;
% %   plot(0:1:48,Q(8,:,4),'-.','LineWidth',5);hold on;
% %   plot(0:1:48,Q(9,:,4),'-','LineWidth',5);hold on;
% %   plot(0:1:48,Q(10,:,4),'--','LineWidth',5);hold on;
% %   plot(0:1:48,Q(11,:,4),':','LineWidth',5);hold on;
% %   title('pS473-AKT','FontWeight','Bold','FontSize',18);
% %   set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %   xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %   ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% %   legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% %   set(gca,'FontWeight','Bold','FontSize',18);
% %   set(gca,'XLim',[0 48]);
% % 
% %  
% %  figure
% %  for i=1:1:11
% %      plot(0:1:48,Q(i,:,1));hold on;
% %  end
% %  title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
% %  set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %  xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %  ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% % legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% % set(gca,'FontWeight','Bold','FontSize',18);
% % figure
% %  for i=1:1:11
% %      plot(0:1:48,Q(i,:,2));title('p-ERK');hold on;
% %  end
% %  title('p-ERK','FontWeight','Bold','FontSize',18);
% %  set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %  set(gca,'FontWeight','Bold','FontSize',18);
% %  xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %  ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% % legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% % set(gca,'FontWeight','Bold','FontSize',18);
% % figure
% %  for i=1:1:11
% %      plot(0:1:48,Q(i,:,3));hold on;
% %  end
% %  ;title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18)
% %  set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %  set(gca,'FontWeight','Bold','FontSize',10);
% %  xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %  ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% % legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% % set(gca,'FontWeight','Bold','FontSize',18);
% % figure
% %  for i=1:1:11
% %      plot(0:1:48,Q(i,:,4));title('pS473-AKT');hold on;
% %  end
% %  title('pS473-AKT','FontWeight','Bold','FontSize',18);
% %  set(gca,'XTick',0:4:48);
% % %  set(gca,'Yticklabel',tick1);
% %  set(gca,'FontWeight','Bold','FontSize',18);
% %  xlabel('Time(hours)','FontWeight','Bold','FontSize',18);
% %  ylabel('Relative value','FontWeight','Bold','FontSize',18);
% % %  set(gca,'XTickLabelRotation',45);
% % legend('Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906');
% % set(gca,'FontWeight','Bold','FontSize',18);

 

