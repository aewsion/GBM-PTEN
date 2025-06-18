close all
A=load('0.1.mat');
B=load('0.5.mat');
C=load('1.mat');
D=load('2.mat');
E=load('0.mat');

label={'0.1','0.5','1','2','Control'};
figure
plot(0:1:200,A.C_T_mean(:),'-','LineWidth',3,'Color','#4682B4');hold on;
plot(0:1:200,B.C_T_mean(:),'-','LineWidth',3,'Color','#778899');hold on;
plot(0:1:200,C.C_T_mean(:),'-','LineWidth',3,'Color','#BC8F8F');hold on;
plot(0:1:200,D.C_T_mean(:),'-','LineWidth',3,'Color','#228B22');hold on;
plot(0:1:200,E.C_T_mean(:),'-','LineWidth',3,'Color','k');hold on;

title('CSF1R-I (WT)','FontWeight','Bold','FontSize',16)
set(gca,'xtick',(0:50:200));
ylim([0 1]);
legend(label,'box','off','Location','Southeast');
xlabel('Days','FontWeight','Bold','FontSize',16);
ylabel('Average Cancer Cell Density','FontWeight','Bold','FontSize',16);
set(gca,'FontWeight','Bold','FontSize',14);
hold on

i=1;
set(i,'Units','Inches');
pos = get(i,'Position');
set(i,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename=['E:\aa文件\project\project2\5\',num2str(5)];
print(i,filename,'-dpdf','-r2000','-r0')