
clc
clear all
close all
global Ptest1
load('pa.mat')
Ptest1=Ptest1;
 Conditions_ODES_2;
     S1=Y1([1:5],[4,5,7]);
    S2=Y2([1,2,3,5,9],[4,5,7]);
    S3=Y1([6:14],[5,6,7]);
    S4=Y3([1,2,3,4,7,13,19,37,73],[5,6,7]);
 sum1=sum(sum((Y4(:,1)-pY1068_EGFR).^2))+sum(sum((Y4(:,2)-p_ERK).^2))+sum(sum((Y4(:,3)-p_IGF1Rbeta).^2))+sum(sum((Y4(:,4)-pS473_AKT).^2))+sum(sum((S1-Exp1).^2))+sum(sum((S2-Exp2).^2))+sum(sum((S3-Exp3).^2))+sum(sum((S4-Exp4).^2))   % y is simulated data; y0 is experimental data
 rmse=sqrt(sum1/84)



    ylabel1={'AKT','GSK3\beta^{Ser9}','Gal-9'};
    ylabel2={'GSK3\beta^{Ser9}','IRF1','Gal-9'};
    tick1={'0 \mum','1 \mum','2 \mum','3 \mum','4 \mum'};
    tick2={'0 h','6 h','12 h','24 h','48 h'};
    tick3={'0 \mum','0.08 \mum','0.16 \mum','0.31 \mum','0.62 \mum','1.25 \mum','2.5 \mum','5 \mum','10 \mum'};
    tick4={'0 min','20 min','40 min','1 h','2 h','4 h','6 h','12 h','24 h'};
    b2=[40/255,120/255,181/255];
    b1=[153/255,77/255,82/255];
    
    for i=1:3
        rmse1=round(sqrt(sum(sum((S1(:,i)-Exp1(:,i)).^2))/5),-4);
        figure,
        plot(1:5,S1(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:5,Exp1(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        axis([1 5 0 1 ]);
        xlabel('Conditions (AKTI 24h)','FontWeight','Bold','FontSize',18)
        xticks(1:1:5)
        set(gca,'Xticklabel',tick1,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
       % ylabel(ylabel1(i),'FontWeight','Bold','FontSize',18);
        title(ylabel1(i),'FontWeight','Bold','FontSize',18);
        legend('Model','Exp');
        txt=['RMSE:',num2str(rmse1)];
        text(4.1,0.8,txt,'FontWeight','Bold')
    end
    
    for i=1:3
        rmse1=round(sqrt(sum(sum((S2(:,i)-Exp2(:,i)).^2))/5),-4);
        figure,
        plot(1:5,S2(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:5,Exp2(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        axis([1 5 0 2.2]);
        xlabel('Conditions (AKTI 2 \mum)','FontWeight','Bold','FontSize',18)
        xticks(1:1:5)
        set(gca,'Xticklabel',tick2,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        %ylabel(ylabel1(i),'FontWeight','Bold','FontSize',18);
        title(ylabel1(i),'FontWeight','Bold','FontSize',18);
        legend('Model','Exp');
        txt=['RMSE: ',num2str(rmse1)];
        text(4.1,0.8,txt,'FontWeight','Bold')
    end
    
    for i=1:3
        rmse1=round(sqrt(sum(sum((S3(:,i)-Exp3(:,i)).^2))/9),-4);
        figure,
        plot(1:9,S3(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:9,Exp3(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        axis([1 9 0 1]);
        xlabel('Conditions (GSK3\betaI 24h)','FontWeight','Bold','FontSize',18)
        xticks(1:1:9)
        set(gca,'Xticklabel',tick3,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        %ylabel(ylabel2(i),'FontWeight','Bold','FontSize',18);
        title(ylabel2(i),'FontWeight','Bold','FontSize',18);
        legend('Model','Exp','location','SouthEast');
        txt=['RMSE: ',num2str(rmse1)];
        text(7,0.2,txt,'FontWeight','Bold')
    end
    
    for i=1:3
        rmse1=round(sqrt(sum(sum((S4(:,i)-Exp4(:,i)).^2))/9),-4);
        figure,
        plot(1:9,S4(:,i),'-o','color',b1,'MarkerFaceColor',b1);hold on;
        plot(1:9,Exp4(:,i),'-o','color',b2,'MarkerFaceColor',b2);
        axis([1 9 0 1]);
        xlabel('Conditions (GSK3\betaI 5 \mum)','FontWeight','Bold','FontSize',18)
        xticks(1:1:9);
        set(gca,'Xticklabel',tick4,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        %ylabel(ylabel2(i),'FontWeight','Bold','FontSize',18);
        title(ylabel2(i),'FontWeight','Bold','FontSize',18);
        legend('Model','Exp','location','SouthEast');
        txt=['RMSE: ',num2str(rmse1)];
        text(7,0.2,txt,'FontWeight','Bold')
    end
    
    tick={'Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906'};
    
    rmse1=round(sqrt(sum(sum((Y4(:,1)-pY1068_EGFR).^2))/11),-4);
    figure,plot(1:11,Y4(:,1),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,pY1068_EGFR,'-o','color',b2,'MarkerFaceColor',b2);
    %set(gca,'ytick',[0.1:2]);
    axis([0 12 0 0.5]);
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    %ylabel('pY1068-EGFR','FontWeight','Bold','FontSize',18);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    %set(gca,'FontWeight','Bold','FontSize',18);
    title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',10);
    legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(5,0.4,txt,'FontWeight','Bold')
    
    rmse1=round(sqrt(sum(sum((Y4(:,2)-p_ERK).^2))/11),-4);
    figure,plot(1:11,Y4(:,2),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,p_ERK,'-o','color',b2,'MarkerFaceColor',b2);
    axis([0 12 0 0.5]);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    %ylabel('p-ERK','FontWeight','Bold','FontSize',18);
    %set(gca,'FontWeight','Bold','FontSize',18);
    title('p-ERK','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',10);
    legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(3,0.4,txt,'FontWeight','Bold')
    
    rmse1=round(sqrt(sum(sum((Y4(:,3)-p_IGF1Rbeta).^2))/11),-4);
    figure,plot(1:11,Y4(:,3),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,p_IGF1Rbeta,'-o','color',b2,'MarkerFaceColor',b2);
    axis([0 12 0 0.5]);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    %ylabel('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
    %set(gca,'FontWeight','Bold','FontSize',18);
    title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',10);
        legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(3,0.4,txt,'FontWeight','Bold')
    
    rmse1=round(sqrt(sum(sum((Y4(:,4)-pS473_AKT).^2))/11),-4);
    figure,plot(1:11,Y4(:,4),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:11,pS473_AKT,'-o','color',b2,'MarkerFaceColor',b2);
    axis([0 12 0 0.5]);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    %ylabel('pS473-AKT','FontWeight','Bold','FontSize',18);
    %set(gca,'FontWeight','Bold','FontSize',18);
    title('pS473-AKT','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',10);
        legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(3,0.4,txt,'FontWeight','Bold')