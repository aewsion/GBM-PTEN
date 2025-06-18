avg_CT=zeros(9,9);
tmp=zeros(37,1);


for i=0:36
    filename=[num2str(i),'.mat'];
    b=load(filename);
    tmp(i+1)=round(mean(b.C_T_mean),2);
end
for i=1:9
    avg_CT(1,i)=tmp(i);
end

a=9;
for i=2:8
    for j=i+1:9
        a=a+1;
        avg_CT(i,j)=tmp(a);
    end
end

avg_CT=avg_CT+avg_CT.';
for i=1:8
    filename=['0',num2str(i),'.mat'];
    b=load(filename);
    avg_CT(i+1,i+1)=round(mean(b.C_T_mean),2);
end
avg_CT(1,1)=tmp(1);
x={'Control','CSF1R-I','Gal9-I','AKT-I','EGFR-I','IGF1R-I','ERK-I','GSK3\beta-I','IRF1-I'};
% h=heatmap(x,x,avg_CT);
% h.FontName='Arial';
%  h.Title='Average cancer cell density (LOSS)';
% h.ColorLimits = [0 1];
% h.Colormap = winter;
% % h.Title.fontsize=12;



SHM=SHeatmap((avg_CT - min(min(avg_CT))) / (max(max(avg_CT))-min(min(avg_CT))),'Format','sq');
SHM=SHM.draw();

clim([0,1])
SHM.setText(); 
% colormap("sky");
ax=gca;
ax.XTickLabel=x;
ax.YTickLabel=x;
title(ax,'Average cancer cell density (WT)','FontWeight','bold','FontSize',18);

ax.FontName="Arial";
ax.FontWeight='bold';
% ax.FontSize=14;
for i=1:size(avg_CT,1)
    for j=1:size(avg_CT,2)
            SHM.setTextMN(i,j,'String',{num2str(avg_CT(i,j))},'FontSize',14,'fontname','arial')

    end
end

i=1;
set(i,'Units','Inches');
pos = get(i,'Position');
set(i,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename=['E:\aa文件\project\project2\5\',num2str(1)];
print(i,filename,'-dpdf','-r2000','-r0')

% dose={'(1,7.4)','(1,50)','(1,9.8)','(1,150)','(25)','(1,160)','(1,11.3)','(1,180)','(1,16)','(500)','(1,9)','(1,140)','(1,30)','(1,400)','(1,16)','(1,200)'};
% a=1;
% for i=2:size(avg_CT,1)
%     for j=[4,6]
%             SHM.setTextMN(i,j,'String',{num2str(avg_CT(i,j)),string(dose(a))},'FontSize',14,'fontname','arial')
%             a=a+1;
% 
%     end
% end
% for i=1:size(avg_CT,1)
%     for j=1:size(avg_CT,2)
%         if avg_CT(i,j)>.9
%             SHM.setTextMN(i,j,'String','**','FontSize',20)
%             SHM.setPatchMN(i,j,'EdgeColor',[1,0,0],'LineWidth',2)
%         end
%     end
% end