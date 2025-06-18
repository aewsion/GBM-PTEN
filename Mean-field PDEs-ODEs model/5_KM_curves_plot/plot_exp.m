clc
close all
addpath('.\MatSurv-master');
%Survival_No=load('No.csv');
Survival_No=load('Nodrug_Dataset.csv');
Survival_D=load('CSF1R-I.csv');

Time_var=200*ones(1,200);
Event_var=ones(1,200);
Group_var = repmat({'CSF1R-I'}, 200, 1);%前100CSF1R-I,后100no treatment。

%Time_var
for i=1:200
    if i<=9
        Time_var(i)=60;
    elseif 9<i&&i<=9+25
        Time_var(i)=80;
    elseif 9+25<i && i<=9+25+21
        Time_var(i)=100;
    elseif 9+25+21<i && i<=9+25+21+1
        Time_var(i)=120;

    elseif 100<i && i<=100+66
        Time_var(i)=20;
    elseif 100+66<i && i<=200
        Time_var(i)=60;
    end

end
%Event_var
for i=56:100%3+29+2：100之间的人存活
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


figure;
hold on
set(gca,'FontSize',15);
plot(Survival_No(:,1),Survival_No(:,2),'color','#34B34A',LineWidth=2.5);hold on;
plot(Survival_D(:,1),Survival_D(:,2),'color','#3855A5',LineWidth=2.5);hold on;
set(gca,'FontSize',14,'FontWeight','bold','fontname','Arial');
xlabel('Time (Days)','FontWeight','Bold','FontSize',16,'fontname','Arial');
ylabel('Percent Survival(%)','FontWeight','Bold','FontSize',16,'fontname','Arial');
title('Experimental Data','FontWeight','Bold','FontSize',16,'fontname','Arial');
legend('No treatment','CSF1R-I','fontname','Arial','box','off');
ylim([0 100])
xlim([0 140])
text(0.63, 0.75, ['p-value = ', num2str(p2 ,'%.3g')], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',12,'FontWeight','bold','fontname','Arial');
box on;
hold off;

i=2;
set(i,'Units','Inches');
pos = get(i,'Position');
set(i,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename=['E:\aa文件\project\project2\4\',num2str(1)];
print(i,filename,'-dpdf','-r2000','-r0')