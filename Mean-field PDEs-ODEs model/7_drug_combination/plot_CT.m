a=load('CT_Loss.mat');
b=load('CT_WT.mat');
CT_WT=a.C_T_mean;
CT_Loss=b.C_T_mean;
tmp1=0.04*19;

x_WT=[10 15 23];
y_WT=[tmp1/19 tmp1/19*1.5 tmp1/19*6];

x_Loss=[8 22];
y_Loss=[tmp1/19 tmp1];
err=[0.0269 0.1704];
% err1=[0.03 0.03*1.5 ];

figure
plot(0:1:30,CT_WT(1:31),'-','LineWidth',3,'Color',[0.85,0.33,0.10]);hold on;
plot(x_WT,y_WT,'^','LineWidth',4,'Color',[0.85,0.33,0.10]);hold on;
plot(0:1:30,CT_Loss(1:31),'-','LineWidth',3,'Color',[0.00,0.45,0.74]);hold on;
% errorbar(x_WT,y_WT,err1,"^","MarkerSize",11,...
    % "MarkerEdgeColor",[0.85,0.33,0.10],"MarkerFaceColor",[0.85,0.33,0.10],'LineWidth',1.5,'Color',[0.85,0.33,0.10])
errorbar(x_Loss,y_Loss,err,"s","MarkerSize",11,...
    "MarkerEdgeColor",[0.00,0.45,0.74],"MarkerFaceColor",[0.00,0.45,0.74],'LineWidth',1.5,'Color',[0.00,0.45,0.74])
% plot(x_Loss,y_Loss,'s','LineWidth',5,'Color',[0.00,0.45,0.74])
% title('PTEN\_WT','FontWeight','Bold','FontSize',18)
legend('Simu. PTEN\_WT','Exp. PTEN\_WT','Simu. PTEN\_Loss','Exp. PTE\_Loss','Location','northwest','fontsize',12)
set(gca,'xtick',(0:5:30));
xlabel('Days','FontWeight','Bold','FontSize',18);
ylabel('Average Cancer Cell Density','FontWeight','Bold','FontSize',18);
set(gca,'FontWeight','Bold','FontSize',16);
hold on