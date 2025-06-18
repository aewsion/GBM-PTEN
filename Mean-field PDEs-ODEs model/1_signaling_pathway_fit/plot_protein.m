
clc
clear all
close all
set(groot, 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial', 'DefaultAxesFontWeight', 'normal', 'DefaultTextFontWeight', 'normal');
global Ptest1
load('pa.mat')
Ptest1=Ptest1;
Conditions_ODES_2;
S1=Y1(1:21,[4,5,7]);%MK2206 24h  剂量[0:0.2:4\muM]
Sim1=Y1([1 1/0.2+1 2/0.2+1 3/0.2+1 4/0.2+1],[4,5,7]);

S2=Y2(:,[4,5,7]);%MK2206 2 \muM 时间[0:1/80:2day]
Sim2=Y2([1,1+20,1+2*20,1+4*20,1+8*20],[4,5,7]);


S3=Y1(22:30,[5,6,7]);%Tideglusib 24h 见31行tick3


S4=Y3(:,[5,6,7]);%Tideglusib 5\muM 时间[0:1/72:1day]
Sim4=Y3([1,2,3,4,7,14,19,37,73],[5,6,7]);


sum1=sum(sum((Y4(:,1)-pY1068_EGFR).^2))+sum(sum((Y4(:,2)-p_ERK).^2))+sum(sum((Y4(:,3)-p_IGF1Rbeta).^2))+sum(sum((Y4(:,4)-pS473_AKT).^2))+sum(sum((Sim1-Exp1).^2))+sum(sum((Sim2-Exp2).^2))+sum(sum((S3-Exp3).^2))+sum(sum((Sim4-Exp4).^2))   % y is simulated data; y0 is experimental data
rmse=sqrt(sum1/84)



ylabel1={'AKT','GSK3\beta^{Ser9}','Gal-9'};
ylabel2={'GSK3\beta^{Ser9}','IRF1','Gal-9'};
tick1={'0 \muM','1 \muM','2 \muM','3 \muM','4 \muM'};
tick2={'0 h','6 h','12 h','24 h','48 h'};
tick3={'0 \muM','0.08 \muM','0.16 \muM','0.31 \muM','0.62 \muM','1.25 \muM','2.5 \muM','5 \muM','10 \muM'};
tick4={'0 min','20 min','40 min','1 h','2 h','4 h','6 h','12 h','24 h'};
b2=[40/255,120/255,181/255];
b1=[153/255,77/255,82/255];

fignum=0;

for i=1:3
    fignum=fignum+1;
    rmse1=round(sqrt(sum(sum((Sim1(:,i)-Exp1(:,i)).^2))/5),4);
    [r,P]=corr(Sim1(:,i),Exp1(:,i));
    figure,
    plot(1:0.2:5,S1(:,i),'-','color','#4169E1','MarkerFaceColor','#4169E1',LineWidth=1);
    hold on;
    plot(1:5,Exp1(:,i),'s','color','#B22222','MarkerFaceColor','#B22222',LineWidth=1);
    axis([1 5 0 1 ]);
    pos=axis;
    
    xticks(1:1:5)
    set(gca,'Xticklabel',tick1,'FontWeight','normal','FontSize',8);
    set(gca,'XTickLabelRotation',45);
    xlabel('Conditions (MK2206 24h)','FontSize',10)
    title(ylabel1(i),'FontSize',12);
    txt=['RMSE:',num2str(rmse1)];
    txt=[txt newline '\rho:' num2str(round(r,4))];
    text(3.8,0.85,txt,'FontSize',10)
    set(gca,'Position',[0.1,0.24,0.8,0.68])

end

tmp1=1:2/40:3;
tmp2=(3+1/40):1/40:4;
tmp3=(4+1/80):1/80:5;
tmp4=[tmp1 tmp2 tmp3];

fignum=fignum+1;
rmse1=round(sqrt(sum(sum((Sim2(:,1)-Exp2(:,1)).^2))/5),4);
[r,P]=corr(Sim2(:,1),Exp2(:,1));
figure,
plot(tmp4,S2(:,1),'-','color','#4169E1','MarkerFaceColor','#4169E1',LineWidth=1);
hold on;
plot(1:5,Exp2(:,1),'s','color','#B22222','MarkerFaceColor','#B22222',LineWidth=1);
axis([1 5 0 2.2]);
xticks(1:1:5)
set(gca,'Xticklabel',tick2,'FontSize',8);
set(gca,'XTickLabelRotation',45);
xlabel('Conditions (MK2206 2 \muM)','FontSize',10)
title(ylabel1(1),'FontSize',12);
leg=legend('Model','Exp','box','off','FontSize',10);
set(leg,'position',[0.64,0.6,0.2,0.1])
txt=['RMSE: ',num2str(rmse1)];
txt=[txt newline '\rho:' num2str(round(r,4))];
text(3.8,0.85*2.2,txt,'FontSize',10)
set(gca,'Position',[0.1,0.24,0.8,0.68])

for i=2:3
    fignum=fignum+1;
    rmse1=round(sqrt(sum(sum((Sim2(:,i)-Exp2(:,i)).^2))/5),4);
    [r,P]=corr(Sim2(:,i),Exp2(:,i));
    figure,
    plot(tmp4,S2(:,i),'-','color','#4169E1','MarkerFaceColor','#4169E1',LineWidth=1);
    hold on;
    plot(1:5,Exp2(:,i),'s','color','#B22222','MarkerFaceColor','#B22222',LineWidth=1);
    axis([1 5 0 1]);
    xticks(1:1:5)
    set(gca,'Xticklabel',tick2,'FontSize',8);
    set(gca,'XTickLabelRotation',45);
    xlabel('Conditions (MK2206 2 \muM)','FontSize',10)
    title(ylabel1(i),'FontSize',12);
    txt=['RMSE: ',num2str(rmse1)];
    txt=[txt newline '\rho:' num2str(round(r,4))];
    text(3.8,0.85,txt,'FontSize',10)
    set(gca,'Position',[0.1,0.24,0.8,0.68])
end


for i=1:3
    fignum=fignum+1;
    rmse1=round(sqrt(sum(sum((S3(:,i)-Exp3(:,i)).^2))/9),4);
    [r,P]=corr(S3(:,i),Exp3(:,i));
    figure,
    plot(1:9,S3(:,i),'-','color','#4169E1','MarkerFaceColor','#4169E1',LineWidth=1);
    hold on;
    plot(1:9,Exp3(:,i),'s','color','#B22222','MarkerFaceColor','#B22222',LineWidth=1);
    axis([1 9 0 1]);
    xticks(1:1:9)
    set(gca,'Xticklabel',tick3,'FontSize',8);
    set(gca,'XTickLabelRotation',45);
    xlabel('Conditions (Tideglusib 24h)','FontSize',10)
    title(ylabel2(i),'FontSize',12);
    txt=['RMSE: ',num2str(rmse1)];
    txt=[txt newline '\rho:' num2str(round(r,4))];
    text(6.5,0.16,txt,'FontSize',10)
    set(gca,'Position',[0.1,0.24,0.8,0.68])
end

%figure中横坐标分布不均匀
tmp1=1:1:4;
tmp2=(4+1/3):1/3:5;
tmp3=(5+1/6):1/6:7;
tmp4=(7+1/18):1/18:8;
tmp5=(8+1/36):1/36:9;
tmp6=[tmp1 tmp2 tmp3 tmp4 tmp5];

for i=1:3
    fignum=fignum+1;
    rmse1=round(sqrt(sum(sum((Sim4(:,i)-Exp4(:,i)).^2))/9),4);
    [r,P]=corr(Sim4(:,i),Exp4(:,i));
    figure,
    plot(tmp6,S4(:,i),'-','color','#4169E1','MarkerFaceColor','#4169E1',LineWidth=1);
    hold on;
    plot(1:9,Exp4(:,i),'s','color','#B22222','MarkerFaceColor','#B22222',LineWidth=1);
    axis([1 9 0 1]);
    xticks(1:1:9);
    set(gca,'Xticklabel',tick4,'FontSize',8);
    set(gca,'XTickLabelRotation',45);
    xlabel('Conditions (Tideglusib 5 \muM)','FontSize',10)
    title(ylabel2(i),'FontSize',12);
    txt=['RMSE: ',num2str(rmse1)];
    txt=[txt newline '\rho:' num2str(round(r,4))];
    text(6.5,0.16,txt,'FontSize',10)
    set(gca,'Position',[0.1,0.24,0.8,0.68])
end

tick={'Control-1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control-2','Gefitinib','Gefitinib 5 \muM','OSI-906','OSI-906 5 \muM'};

fignum=fignum+1;
rmse1=round(sqrt(sum(sum((Y4(:,1)-pY1068_EGFR).^2))/11),4);
[r,p]=corr(Y4(:,1),pY1068_EGFR);
figure,
pBar = bar(1:11, [Y4(:,1), pY1068_EGFR], 'grouped'); 
pBar(1).FaceColor = '#4169E1';
pBar(2).FaceColor = '#B22222';
axis([0 12 0 0.5]);
set(gca,'XTick',1:1:11);
set(gca,'XTickLabelRotation',45);
set(gca,'Xticklabel',tick,'FontSize',8);
xlabel('Conditions','FontSize',10);
title('pY1068-EGFR','FontSize',12);
leg=legend('Model','Exp','box','off','FontSize',10);
set(leg,'position',[0.1,0.8,0.2,0.1])
txt=['RMSE: ',num2str(rmse1)];
txt=[txt newline '\rho:' num2str(round(r,4))];
text(8.5,0.43,txt,'FontSize',10)
set(gca,'Position',[0.1,0.33,0.8,0.6])

fignum=fignum+1;
rmse1=round(sqrt(sum(sum((Y4(:,2)-p_ERK).^2))/11),4);
[r,p]=corr(Y4(:,2),p_ERK);
figure,
pBar = bar(1:11, [Y4(:,2), p_ERK], 'grouped'); 
pBar(1).FaceColor = '#4169E1';
pBar(2).FaceColor = '#B22222';
axis([0 12 0 0.5]);
set(gca,'XTick',1:1:11);
set(gca,'XTickLabelRotation',45);
set(gca,'Xticklabel',tick,'FontSize',8);
xlabel('Conditions','FontSize',10);
title('p-ERK','FontSize',12);
txt=['RMSE: ',num2str(rmse1)];
txt=[txt newline '\rho:' num2str(round(r,4))];
text(8.5,0.43,txt,'FontSize',10)
set(gca,'Position',[0.1,0.33,0.8,0.6])

fignum=fignum+1;
rmse1=round(sqrt(sum(sum((Y4(:,3)-p_IGF1Rbeta).^2))/11),4);
[r,p]=corr(Y4(:,3),p_IGF1Rbeta);
figure,
pBar = bar(1:11, [Y4(:,3), p_IGF1Rbeta], 'grouped'); 
pBar(1).FaceColor = '#4169E1';
pBar(2).FaceColor = '#B22222';
axis([0 12 0 0.5]);
set(gca,'XTick',1:1:11);
set(gca,'XTickLabelRotation',45);
set(gca,'Xticklabel',tick,'FontSize',8);
xlabel('Conditions','FontSize',10);
title('p-IGF1R\beta','FontSize',12);
txt=['RMSE: ',num2str(rmse1)];
txt=[txt newline '\rho:' num2str(round(r,4))];
text(8.5,0.43,txt,'FontSize',10)
set(gca,'Position',[0.1,0.33,0.8,0.6])

fignum=fignum+1;
rmse1=round(sqrt(sum(sum((Y4(:,4)-pS473_AKT).^2))/11),4);
[r,p]=corr(Y4(:,4),pS473_AKT);
figure,
pBar = bar(1:11, [Y4(:,4), pS473_AKT], 'grouped'); 
pBar(1).FaceColor = '#4169E1';
pBar(2).FaceColor = '#B22222';
axis([0 12 0 0.5]);
set(gca,'XTick',1:1:11);
set(gca,'XTickLabelRotation',45);
set(gca,'Xticklabel',tick,'FontSize',8);
xlabel('Conditions','FontSize',10);
title('pS473-AKT','FontSize',12);
txt=['RMSE: ',num2str(rmse1)];
txt=[txt newline '\rho:' num2str(round(r,4))];
text(8.5,0.43,txt,'FontSize',10)
set(gca,'Position',[0.1,0.33,0.8,0.6])

figHandles = findobj('Type', 'figure');
for i = 1:length(figHandles)
    fig = figHandles(i);
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    filename=strcat("C:\Users\Administrator\OneDrive\Desktop\YangH\code_figures\上传\figures_20250209\",num2str(i),".pdf");
    print(fig, filename,'-dpdf','-vector', '-r300')
end
