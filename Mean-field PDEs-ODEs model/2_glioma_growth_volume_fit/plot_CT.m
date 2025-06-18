clc
close all
a=load('CT_Loss.mat');
b=load('CT_WT.mat');
CT_WT=a.C_T_mean;
CT_Loss=b.C_T_mean;
tmp1=0.04*19;

x_WT=[10 15 23];
y_WT=[tmp1/19 tmp1/19*1.5 tmp1/19*6];

x_Loss=[8 22];
y_Loss=[tmp1/19 tmp1/19*19.6];
err=[0.0269 0.1704];
% err1=[0.03 0.03*1.5 ];

figure
plot(0:1:30,CT_WT(1:31),'-','LineWidth',3,'Color',[0.85,0.33,0.10]);hold on;
plot(x_WT,y_WT,'^','LineWidth',4,'Color',[0.85,0.33,0.10]);hold on;
legend('Simu. PTEN\_WT','Exp. PTEN\_WT','Location','northwest','fontsize',12,'box','off')
set(gca,'xtick',(0:5:30));
xlabel('Days','FontWeight','Bold','FontSize',18);
ylabel('Average Cancer Cell Density','FontWeight','Bold','FontSize',18);
set(gca,'FontWeight','Bold','FontSize',16);
ylim([0 1])
hold on

figure
plot(0:1:30,CT_Loss(1:31),'-','LineWidth',3,'Color',[0.00,0.45,0.74]);hold on;
errorbar(x_Loss,y_Loss,err,"s","MarkerSize",11,...
    "MarkerEdgeColor",[0.00,0.45,0.74],"MarkerFaceColor",[0.00,0.45,0.74],'LineWidth',1.5,'Color',[0.00,0.45,0.74])
legend('Simu. PTEN\_Null','Exp. PTE\_Null','Location','northwest','fontsize',12,'box','off')
set(gca,'xtick',(0:5:30));
xlabel('Days','FontWeight','Bold','FontSize',18);
ylabel('Average Cancer Cell Density','FontWeight','Bold','FontSize',18);
set(gca,'FontWeight','Bold','FontSize',16);
hold on

% for i=1:2
%     set(i,'Units','Inches');
%     pos = get(i,'Position');
%     set(i,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     filename=['E:\aa文件\project\project2\2\',num2str(i)];
%     print(i,filename,'-dpdf','-r2000','-r0')
% end