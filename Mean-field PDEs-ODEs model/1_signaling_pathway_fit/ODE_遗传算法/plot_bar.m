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
        rmse1=roundn(sqrt(sum(sum((S1(:,i)-Exp1(:,i)).^2))/5),-4);
        y=[S1(1,i),Exp1(1,i);S1(2,i),Exp1(2,i);S1(3,i),Exp1(3,i);S1(4,i),Exp1(4,i);S1(5,i),Exp1(5,i)];
        figure,
        h=bar(1:5,y);
        axis([0 6 0 1 ]);
        set(h(1),'FaceColor',b1); 
        set(h(2),'FaceColor',b2); 
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:5)
        set(gca,'Xticklabel',tick1,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        ylabel(ylabel1(i),'FontWeight','Bold','FontSize',18);
        title('AKTI 24h','FontWeight','Bold','FontSize',18);
        legend('Model','Exp');
        txt=['RMSE:',num2str(rmse1)];
        text(4.7,0.8,txt,'FontWeight','Bold')
    end
    
    for i=1:3
        rmse1=roundn(sqrt(sum(sum((S2(:,i)-Exp2(:,i)).^2))/5),-4);
        y=[S2(1,i),Exp2(1,i);S2(2,i),Exp2(2,i);S2(3,i),Exp2(3,i);S2(4,i),Exp2(4,i);S2(5,i),Exp2(5,i)];
        figure,
        h=bar(1:5,y);
        axis([0 6 0 1 ]);
        set(h(1),'FaceColor',b1); 
        set(h(2),'FaceColor',b2); 
        axis([0 6 0 1]);
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:5)
        set(gca,'Xticklabel',tick2,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        ylabel(ylabel1(i),'FontWeight','Bold','FontSize',18);
        title('AKTI 2 \mum','FontWeight','Bold','FontSize',18);
        legend('Model','Exp');
        txt=['RMSE: ',num2str(rmse1)];
        text(4.7,0.8,txt,'FontWeight','Bold')
    end
    
    for i=1:3
        rmse1=roundn(sqrt(sum(sum((S3(:,i)-Exp3(:,i)).^2))/5),-4);
         y=[S3(1,i),Exp3(1,i);S3(2,i),Exp3(2,i);S3(3,i),Exp3(3,i);S3(4,i),Exp3(4,i);S3(5,i),Exp3(5,i);S3(6,i),Exp3(6,i);S3(7,i),Exp3(7,i);S3(8,i),Exp3(8,i);S3(9,i),Exp3(9,i)];
        figure,
        h=bar(1:9,y);
        axis([0 10 0 1 ]);
        set(h(1),'FaceColor',b1); 
        set(h(2),'FaceColor',b2); 
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:9)
        set(gca,'Xticklabel',tick3,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        ylabel(ylabel2(i),'FontWeight','Bold','FontSize',18);
        title('GSK3\betaI 24h','FontWeight','Bold','FontSize',18);
        legend('Model','Exp','location','NorthWest');
        txt=['RMSE: ',num2str(rmse1)];
        text(0.2,0.8,txt,'FontWeight','Bold')
    end
    
    for i=1:3
        rmse1=roundn(sqrt(sum(sum((S3(:,i)-Exp3(:,i)).^2))/5),-4);
        y=[S4(1,i),Exp4(1,i);S4(2,i),Exp4(2,i);S4(3,i),Exp4(3,i);S4(4,i),Exp4(4,i);S4(5,i),Exp4(5,i);S4(6,i),Exp4(6,i);S4(7,i),Exp4(7,i);S4(8,i),Exp4(8,i);S4(9,i),Exp4(9,i)];
        figure,
        h=bar(1:9,y);
        axis([0 10 0 1 ]);
        set(h(1),'FaceColor',b1); 
        set(h(2),'FaceColor',b2); 
        xlabel('Conditions','FontWeight','Bold','FontSize',18)
        xticks(1:1:9);
        set(gca,'Xticklabel',tick4,'FontWeight','Bold','FontSize',10);
        set(gca,'XTickLabelRotation',45);
        ylabel(ylabel2(i),'FontWeight','Bold','FontSize',18);
        title('GSK3\betaI 5 \mum','FontWeight','Bold','FontSize',18);
        legend('Model','Exp','location','NorthWest');
        txt=['RMSE: ',num2str(rmse1)];
        text(0.2,0.8,txt,'FontWeight','Bold')
    end
    
    tick={'Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906'};
    
    rmse1=roundn(sqrt(sum(sum((Y4(:,1)-pY1068_EGFR).^2))/11),-4);
    y=[Y4(1,1),pY1068_EGFR(1);Y4(2,1),pY1068_EGFR(2);Y4(3,1),pY1068_EGFR(3);Y4(4,1),pY1068_EGFR(5);Y4(5,1),pY1068_EGFR(5);Y4(6,1),pY1068_EGFR(6);Y4(7,1),pY1068_EGFR(7);Y4(8,1),pY1068_EGFR(8);Y4(9,1),pY1068_EGFR(9);Y4(10,1),pY1068_EGFR(10);Y4(11,1),pY1068_EGFR(11)];
    figure,
    h=bar(1:11,y);
    axis([0 13 0 0.5 ]);
    set(h(1),'FaceColor',b1); 
    set(h(2),'FaceColor',b2); 
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    ylabel('pY1068-EGFR','FontWeight','Bold','FontSize',18);
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    set(gca,'FontWeight','Bold','FontSize',18);
%     title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
    legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(5,0.4,txt,'FontWeight','Bold')
    
    rmse1=roundn(sqrt(sum(sum((Y4(:,2)-p_ERK).^2))/11),-4);
    y=[Y4(1,2),p_ERK(1);Y4(2,2),p_ERK(2);Y4(3,2),p_ERK(3);Y4(4,2),p_ERK(5);Y4(5,2),p_ERK(5);Y4(6,2),p_ERK(6);Y4(7,2),p_ERK(7);Y4(8,2),p_ERK(8);Y4(9,2),p_ERK(9);Y4(10,2),p_ERK(10);Y4(11,2),p_ERK(11)];
    figure,
    h=bar(1:11,y);
    axis([0 13 0 0.5 ]);
    set(h(1),'FaceColor',b1); 
    set(h(2),'FaceColor',b2); 
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    ylabel('p-ERK','FontWeight','Bold','FontSize',18);
    set(gca,'FontWeight','Bold','FontSize',18);
%     title('p-ERK','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
    legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(3,0.4,txt,'FontWeight','Bold')
    
    rmse1=roundn(sqrt(sum(sum((Y4(:,3)-p_IGF1Rbeta).^2))/11),-4);
    y=[Y4(1,3),p_IGF1Rbeta(1);Y4(2,3),p_IGF1Rbeta(2);Y4(3,3),p_IGF1Rbeta(3);Y4(4,3),p_IGF1Rbeta(5);Y4(5,3),p_IGF1Rbeta(5);Y4(6,3),p_IGF1Rbeta(6);Y4(7,3),p_IGF1Rbeta(7);Y4(8,3),p_IGF1Rbeta(8);Y4(9,3),p_IGF1Rbeta(9);Y4(10,3),p_IGF1Rbeta(10);Y4(11,3),p_IGF1Rbeta(11)];
    figure,
    h=bar(1:11,y);
    axis([0 13 0 0.5 ]);
    set(h(1),'FaceColor',b1); 
    set(h(2),'FaceColor',b2); 
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    ylabel('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
    set(gca,'FontWeight','Bold','FontSize',18);
%     title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
    legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(3,0.4,txt,'FontWeight','Bold')
    
    rmse1=roundn(sqrt(sum(sum((Y4(:,4)-pS473_AKT).^2))/11),-4);
    y=[Y4(1,4),pS473_AKT(1);Y4(2,4),pS473_AKT(2);Y4(3,4),pS473_AKT(3);Y4(4,4),pS473_AKT(5);Y4(5,4),pS473_AKT(5);Y4(6,4),pS473_AKT(6);Y4(7,4),pS473_AKT(7);Y4(8,4),pS473_AKT(8);Y4(9,4),pS473_AKT(9);Y4(10,4),pS473_AKT(10);Y4(11,4),pS473_AKT(11)];
    figure,
    h=bar(1:11,y);
    axis([0 13 0 0.5 ]);
    set(h(1),'FaceColor',b1); 
    set(h(2),'FaceColor',b2); 
    set(gca,'XTick',1:1:11);
    set(gca,'XTickLabelRotation',45);
    xlabel('Condition','FontWeight','Bold','FontSize',18);
    ylabel('pS473-AKT','FontWeight','Bold','FontSize',18);
    set(gca,'FontWeight','Bold','FontSize',18);
%     title('pS473-AKT','FontWeight','Bold','FontSize',18);
    set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
        legend('Model','Exp');
    txt=['RMSE: ',num2str(rmse1)];
    text(3,0.4,txt,'FontWeight','Bold')