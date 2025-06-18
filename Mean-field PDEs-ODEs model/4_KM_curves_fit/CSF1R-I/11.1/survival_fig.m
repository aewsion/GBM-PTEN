global PTEN
PTEN=1;


if PTEN==0
    load('stochastic_KO.mat');
else
    load('stochastic_WT.mat');
end
[asizeC,bsizeC]=size(C);
TimeLength=bsizeC-1;
threshold=0.951;
[p,p12]=size(C);

%death_day表示死亡的天数
death_day=zeros(p,1);
for person=1:1:p
    if isempty(find(C(person,:)>=threshold,1))%为空进入
        max_death_day=max(C(person,:));
        i=find(C(person,:)==max_death_day);
        if i+50<=p12 && max_death_day>=0.94 && (max_death_day-sum(C(person,i+1:i+50))/50)/max_death_day<0.005
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
%death_person_10表示每10天死亡的人数
death_person_10=zeros(TimeLength/10+1,1);
for i=1:1:TimeLength/10
    death_person_10(i+1,1)=sum(death_person((i-1)*10+1:i*10,1));
end
if PTEN==0
    save('stochastic_KO_death.mat','Survival_Rate','death_person','death_person_10','death_day');
else
    save('stochastic_MT_death.mat','Survival_Rate','death_person','death_person_10','death_day');
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
    xlabel('Time (Days)','FontWeight','Bold','FontSize',18);
    ylabel('Cancer Cell Density','FontWeight','Bold','FontSize',18);
    set(gca,'FontSize',18);
    figure_FontSize=18;
    ylim([0.5 1]);
    title('PTEN-null','FontWeight','Bold','FontSize',18);
    set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    set(findobj('FontSize',18),'FontSize',figure_FontSize);
    box on;

    ST2=[0 60 30 10 0 0 0 0];%20day-90day
    PST2=death_person_10(3:10,1);
    figure,
    plot(20:10:90,ST2,'b o-- ','Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','b');axis([0,140,0,100])
    hold on, 
    plot(20:10:90,PST2,'r o-- ','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r');axis([0,140,0,100])
    set(gca,'xtick',(20:10:90));
    xlim([20 90]);
    legend('Experimental data','Prediction')
    title('PTEN-null','FontWeight','Bold','FontSize',18);
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
        plot(20:90,C(person,21:91));
        % plot(0:90,C(person,1:91));
    end
    set(gca,'xtick',(1:10:90));
    xlim([20 90]);
    xlabel('Time (Days)','FontWeight','Bold','FontSize',18);
    ylabel('Cancer Cell Density','FontWeight','Bold','FontSize',18);
    set(gca,'FontSize',18);
    figure_FontSize=18;
    ylim([0.5 1]);
    set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    % set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    set(findobj('FontSize',18),'FontSize',figure_FontSize);
    title('PTEN-WT','FontWeight','Bold','FontSize',18);
    box on;

    ST1=[0 0 24 36 27 10 3 0];%20day-90day
    PST1=death_person_10(3:10,1);
    figure,
    plot(20:10:90,ST1,'b o-- ','Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','b');axis([0,140,0,100])
    hold on, 
    plot(20:10:90,PST1,'r o-- ','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r');axis([0,140,0,100])
    set(gca,'xtick',(20:10:90));
    xlim([20 90]);
    legend('Experimental data','Prediction')
    title('PTEN-WT','FontWeight','Bold','FontSize',18);
    set(gca,'xticklabel',{'20','30','40','50','60','70','80','90'})
    xlabel('Time (Days)','FontWeight','Bold','FontSize',18);
    ylabel('Frequency','FontWeight','Bold','FontSize',18);
    set(gca,'FontSize',18);
    figure_FontSize=18;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    set(findobj('FontSize',18),'FontSize',figure_FontSize);
end