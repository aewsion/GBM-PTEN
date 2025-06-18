clc
close all
addpath('.\MatSurv-master');
a=load('stochastic_CSF1R-I.mat');
b=load('stochastic_No.mat');

Time_var=200*ones(1,200);
Event_var=ones(1,200);
Group_var = repmat({'CSF1R-I'}, 200, 1);%前100CSF1R-I,后100no treatment。

%Time_var
for i=1:200
    if i<=3
        Time_var(i)=60;
    elseif 3<i&&i<=3+29
        Time_var(i)=80;
    elseif 3+29<i && i<=3+29+20
        Time_var(i)=100;
    elseif 3+29+20<i && i<=3+29+20+2
        Time_var(i)=120;

    elseif 100<i && i<=65+100
        Time_var(i)=20;
    elseif 100+65<i&&i<=100+65+13
        Time_var(i)=60;
    end

end
%Event_var
for i=3+29+20+2:100%3+29+2：100之间的人存活
    Event_var(i)=0;
end
for i=100+65+13:200
    Event_var(i)=0;
end
%Group_var
for i=101:1:200
    Group_var(i) = {'No treatment'};
end
%fig
Time_var = Time_var';
Event_var = Event_var';
[p2,~,~] = MatSurv(Time_var, Event_var, Group_var);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Survival_No=100*ones(1,200);
Survival_D=100*ones(1,200);
death_day=zeros(1,100);
death_person=zeros(1,200);
threshold=0.7;
for person=1:1:100
    if isempty(find(b.C(person,:)>=threshold,1))%为空进入
        max_death_day=max(b.C(person,:));
        i=find(b.C(person,:)==max_death_day);
        if i+50<=201 && max_death_day>=0.77 && (max_death_day-sum(b.C(person,i+1:i+50))/50)/max_death_day<0.005
            aaaaaa=1;
            death_day(person,1)=i;
        end
    else
        death_day(person,1)=find(b.C(person,:)>=threshold,1); %只要第一个，大于阈值的第一天
    end
end


for i=1:1:200+1
    dp=find(death_day==i);
    if isempty(dp)%空就是这天无人死亡
        continue;
    else
        death_person(i,1)=length(dp);
    end
end

for t=1:200
    Survival_No(t)=100-sum(death_person(1:t));
    Survival_D(t)=100-sum(a.C(:,t)>=0.7);
end

figure;
hold on
set(gca,'FontSize',15);
plot(1:140,Survival_No(1:140),'color','#34B34A',LineWidth=2.5);hold on;
plot(1:140,Survival_D(1:140),'color','#3855A5',LineWidth=2.5);hold on;
set(gca,'FontSize',14,'FontWeight','bold','fontname','Arial');
xlabel('Time (Days)','FontWeight','Bold','FontSize',16,'fontname','Arial');
ylabel('Percent Survival(%)','FontWeight','Bold','FontSize',16,'fontname','Arial');
title('Prediction','FontWeight','Bold','FontSize',16,'fontname','Arial');
% legend('No treatment','CSF1R-I','fontname','Arial','box','off');
ylim([0 100])
xlim([0 140])
text(0.63, 0.75, ['p-value = ', num2str(p2 ,'%.3g')], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',12,'FontWeight','bold','fontname','Arial');
box on;
hold off;

i=2;
set(i,'Units','Inches');
pos = get(i,'Position');
set(i,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename=['E:\aa文件\project\project2\4\',num2str(2)];
print(i,filename,'-dpdf','-r2000','-r0')
