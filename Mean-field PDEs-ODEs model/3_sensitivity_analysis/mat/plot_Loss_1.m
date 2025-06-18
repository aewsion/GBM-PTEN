
clear all
close all
load('Loss_1.mat');
PCC_RG=PCC_RG([1 3 18 4 5 6:9 13 15]);
Para_label={'{\alpha}_{AKT}','{\beta}_{C}','{\beta}_{G}','{\beta}_{CI}','{\beta}_{M2M1}','{K_1}^*','{K_2}^*','K_3','K_4','d_{TM1}','{S_A}^*'};
b1=[0,0,0.5451];
num_para=11;


figure,
hh1=bar(1:num_para,PCC_RG);  hold on;
for i=[1 4 11]%[1 2 3]
    text(i,PCC_RG(i)-0.07,'*','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',18);
end
for i=[5 10]%[3 5 10]%[5 6 7 8 10]
    text(i,PCC_RG(i)-0.1,'*','HorizontalAlignment','center','FontSize',18);
end
text(8,0.8,'* P-value <0.01','FontSize',14 ,'FontName','Arial','FontWeight','bold');
set(hh1,'facecolor', b1);
ylim([-1 1])
set(gca,'xtick',1:1:num_para)
set(gca,'xticklabel',Para_label,'FontSize',14,'FontWeight','Bold','FontName','Arial');
set(gca,'XTickLabelRotation',60);
ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',16);
title('PTEN\_Null CSF1R-I','FontWeight','Bold','FontSize',16);
fontname("Arial")

% i=1;
% set(i,'Units','Inches');
% pos = get(i,'Position');
% set(i,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% filename=['E:\aa文件\project\project2\3\',num2str(4)];
% print(i,filename,'-dpdf','-r2000','-r0')