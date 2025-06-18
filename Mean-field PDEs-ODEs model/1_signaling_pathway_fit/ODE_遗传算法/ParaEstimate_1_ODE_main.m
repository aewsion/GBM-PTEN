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

group=1;%想要用遗传算法来估计10组参数
NVAR =46;%model 2            % Number of variables, according to the original setting of Para0  每组参数一共有22个参数
p=ones(group,NVAR+1);     %p来保存10组参数的22个参数加上每一组参数的评估值（计算得到），一共10×23大小
for q=1:1:group
    NIND = 20;           % Number of individuals per subpopulations    相当于每次20条染色体
    MAXGEN = 6500;        % maximum Number of generations     每找一组参数迭代的次数 400次
    GGAP = 0.90;         % Generation gap, how many new individuals are created
    PRECI = 8;           % Precision of binary representation  精确到8位二进制
    Para0= ones(1,NVAR);  % initial values in parameter space
    Ptest1=Para0;
    shu_min=0.00001;           %Minimum search range    首先设置每个参数的最小值
    shu_max=50;          %Maximum search range        再设置有的参数的最大值
    % Build field descriptor
    % FieldD = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,20) 1;shu_max*ones(1,10) 1*ones(1,8) 10*ones(1,2) 11]
    %               rep([1; 0; 1 ;1], [1, NVAR])];%model 1
    FieldD = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,23) 0.25 shu_min*ones(1,2) 0.00005*ones(1,8) ones(1,8) 4 shu_min*ones(1,2) shu_min;shu_max*ones(1,17) 3*ones(1,14) 10*ones(1,3) 10*ones(1,8) 4 3*ones(1,3)]
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
    [value ind]=min(ObjV);
    Ptest1=chrom_value(ind,:)
    Conditions_ODES_2;
    S1=Y1([1:5],[4,5,7]);
    S2=Y2([1,2,3,5,9],[4,5,7]);
    S3=Y1([6:14],[5,6,7]);
    S4=Y3([1,2,3,4,7,13,19,37,73],[5,6,7]);
    
    title1={'AKT AKTI 24h','GSK3\beta AKTI 24h','Gal9 AKTI 24h'};
    title2={'AKT AKTI 2 \mum','GSK3\beta AKTI 2 \mum','Gal9 AKTI 2 \mum'};
    title3={'GSK3\beta GSK3\betaI 24h','IRF1 GSK3\betaI 24h','Gal9 GSK3\betaI 24h'};
    title4={'GSK3\beta GSK3\betaI 5 \mum','IRF1 GSK3\betaI 5 \mum','Gal9 GSK3\betaI 5 \mum'};
    tick1={'0 \mum','1 \mum','2 \mum','3 \mum','4 \mum'};
    tick2={'0 h','6 h','12 h','24 h','48 h'};
    tick3={'0 \mum','0.08 \mum','0.16 \mum','0.31 \mum','0.62 \mum','1.25 \mum','2.5 \mum','5 \mum','10 \mum'};
    tick4={'0 min','20 min','40 min','1 h','2 h','4 h','6 h','12 h','24 h'};
    b1=[147/255,224/255,255/255];
    b2=[153/255,77/255,82/255];
    
    for i=1:3
        figure,
        plot(1:5,S1(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:5,Exp1(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:5)
        set(gca,'Xticklabel',tick1,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        title(title1(i));
        legend('Model','Exp');
    end
    
    for i=1:3
        figure,
        plot(1:5,S2(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:5,Exp2(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:5)
        set(gca,'Xticklabel',tick2,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        title(title2(i));
        legend('Model','Exp');
    end
    
    for i=1:3
        figure,
        plot(1:9,S3(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:9,Exp3(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:9)
        set(gca,'Xticklabel',tick3,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        title(title3(i));
        legend('Model','Exp');
    end
    
    for i=1:3
        figure,
        plot(1:9,S4(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:9,Exp4(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:9);
        set(gca,'Xticklabel',tick4,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        title(title4(i));
        legend('Model','Exp');
    end
    
    tick={'Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906'};
    figure,plot(1:11,Y4(:,1),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,pY1068_EGFR,'-o','color',b2,'MarkerFaceColor',b2);
    %set(gca,'ytick',[0.1:2]);
    %  xlabel('Condition','FontWeight','Bold','FontSize',18);
    ylabel('Relative value','FontWeight','Bold','FontSize',18);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    legend('Model','Experiments');
    set(gca,'FontWeight','Bold','FontSize',18);
    title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
    figure,plot(1:11,Y4(:,2),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,p_ERK,'-o','color',b2,'MarkerFaceColor',b2);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    set(gca,'YTick',0.4:0.1:1);
    set(gca,'Yticklabel',0.4:0.1:1,'FontWeight','Bold','FontSize',18);
    ylabel('Relative value','FontWeight','Bold','FontSize',18);
    legend('Model','Experiments');
    set(gca,'FontWeight','Bold','FontSize',18);
    title('p-ERK','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
    figure,plot(1:11,Y4(:,3),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,p_IGF1Rbeta,'-o','color',b2,'MarkerFaceColor',b2);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    %  xlabel('Condition','FontWeight','Bold','FontSize',18);
    ylabel('Relative value','FontWeight','Bold','FontSize',18);
    legend('Model','Experiments');
    set(gca,'FontWeight','Bold','FontSize',18);
    title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
    figure,plot(1:11,Y4(:,4),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,pS473_AKT,'-o','color',b2,'MarkerFaceColor',b2);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    %  xlabel('Condition','FontWeight','Bold','FontSize',18);
    ylabel('Relative value','FontWeight','Bold','FontSize',18);
    legend('Model','Experiments');
    set(gca,'FontWeight','Bold','FontSize',18);
    title('pS473-AKT','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
end
save('pa.mat','Ptest1');





