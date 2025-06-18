clc
close all
global PTEN
PTEN=1;


if PTEN==0
    load('stochastic_KO.mat');
else
    load('stochastic_WT.mat');
end
[asizeC,bsizeC]=size(C);
TimeLength=bsizeC-1;
threshold=0.7;
[p,p12]=size(C);

%death_day表示死亡的天数
death_day=zeros(p,1);
for person=1:1:p
    if isempty(find(C(person,:)>=threshold,1))%为空进入
        max_death_day=max(C(person,:));
        i=find(C(person,:)==max_death_day);
        if i+50<=p12 && max_death_day>=0.97 && (max_death_day-sum(C(person,i+1:i+50))/50)/max_death_day<0.005
            aaaaaa=1
            death_day(person,1)=i;
        end
    else
        death_day(person,1)=find(C(person,:)>=threshold,1); %只要第一个，大于阈值的第一天
    end
end
%death_person表示每天死亡的人数
death_person=zeros(TimeLength+1,1);
for i=1:1:TimeLength+1
    dp=find(death_day==i);
    if isempty(dp)%空就是这天无人死亡
        continue;
    else
        death_person(i,1)=length(dp);
    end
end
Survival_Rate=100*ones(TimeLength+1,1);
for i=2:1:TimeLength+1
    Survival_Rate(i,1)=(p-death_person(i,1))*100/p;
end
%death_person_20表示每20天死亡的人数
death_person_20=zeros(TimeLength/20+1,1);
for i=1:1:TimeLength/20
    death_person_20(i+1,1)=sum(death_person((i-1)*20+1:i*20,1));
end
if PTEN==1
    save('stochastic_CSF1RI_MT.mat','Survival_Rate','death_person','death_person_20','death_day');
else
    save('stochastic_Gal9I_KO.mat','Survival_Rate','death_person','death_person_10','death_day');
end


if PTEN==0
    %without treatment KO
    figure,
    for person=1:p
        hold on;
         plot(20:90,C(person,21:91));
         % plot(0:90,C(person,1:91));
    end
    set(gca,'xtick',(1:10:90));
    xlim([20 90]);
    set(gca,'FontSize',14,'FontWeight','bold');
    xlabel('Time (Days)','FontWeight','Bold','FontSize',16);
    ylabel('Cancer Cell Density','FontWeight','Bold','FontSize',18);
    % figure_FontSize=18;
    ylim([0.5 1]);
    title('No Treatment','FontWeight','Bold','FontSize',16);
    % set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    % set(findobj('FontSize',18),'FontSize',figure_FontSize);
    box on;

    ST2=[0 60 30 10 0 0 0 0];%20day-90day
    PST2=death_person_10(3:10,1);
    figure,
    plot(0:20:140,ST2,'b o-- ','Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','b');axis([0,140,0,100])
    hold on, 
    plot(0:20:140,PST2,'r o-- ','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r');axis([0,140,0,100])
    set(gca,'xtick',(0:20:140));
    xlim([0 140]);
    legend('Experimental data','Prediction')
    title('PTEN-null CSF1RI','FontWeight','Bold','FontSize',18);
    set(gca,'xticklabel',{'20','30','40','50','60','70','80','90'})
    xlabel('Time (Days)','FontWeight','Bold','FontSize',18);
    ylabel('Frequency','FontWeight','Bold','FontSize',18);
    set(gca,'FontSize',18);
    figure_FontSize=18;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    set(findobj('FontSize',18),'FontSize',figure_FontSize);
else
    
    figure,
    for person=1:p
        hold on;
        plot(0:1:200,C(person,:));
        % plot(0:90,C(person,1:91));
    end
    set(gca,'xtick',(0:50:200));
    xlim([0 200]);
    set(gca,'FontSize',14,'FontWeight','bold','fontname','Arial');
    xlabel('Time (Days)','FontWeight','Bold','FontSize',16,'fontname','Arial');
    ylabel('Cancer Cell Density','FontWeight','Bold','FontSize',16,'fontname','Arial');
    ylim([0 1]);
    title('No treatment','FontWeight','Bold','FontSize',16,'fontname','Arial');
    box on;

    ST1=[0 66 34 0 0 0 0 0];%20day-90day
    PST1=death_person_20(1:8,1);
    figure,
    plot(0:20:140,ST1,'o-- ','Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','b');axis([0,140,0,100])
    hold on, 
    plot(0:20:140,PST1,'o-- ','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r');axis([0,140,0,100])
    % set(gca,'xtick',(0:20:140));
    xlim([0 140]);
    set(gca,'FontSize',14,'FontWeight','bold','fontname','Arial');
    legend('Experimental Data','Prediction','fontname','Arial','box','off')
    title('No treatment','FontWeight','Bold','FontSize',16,'fontname','Arial');
    % set(gca,'xticklabel',{'20','30','40','50','60','70','80','90'})
    xlabel('Time (Days)','FontWeight','Bold','FontSize',16,'fontname','Arial');
    ylabel('Frequency','FontWeight','Bold','FontSize',16,'fontname','Arial');
    % set(gca,'FontSize',18);
    % figure_FontSize=18;
    % set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    % set(findobj('FontSize',18),'FontSize',figure_FontSize);
end

% i=2;
% set(i,'Units','Inches');
% pos = get(i,'Position');
% set(i,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% filename=['E:\aa文件\project\project2\4\',num2str(5)];
% print(i,filename,'-dpdf','-r2000','-r0')