clc
close all
addpath('MatSurv-master');

simulation_no_treatment=load('stochastic_No.mat');
exp_no_treatment=load('Nodrug_Dataset.csv');

Time_var=60*ones(1,130);
Event_var=zeros(1,130);
Group_var = repmat({'prediction'}, 130, 1); % 前100 prediction,后30 experience。
%Group_var
for i=101:1:130
    Group_var(i) = {'experience'};
end

ll = length(exp_no_treatment);
tmp = 1;
for ii=1:1:(ll-1)
    k = exp_no_treatment(ii,3) - exp_no_treatment(ii+1,3);
    if(k > 0)
        for jj = tmp:1:(tmp+k-1)
            Event_var(jj+100) = 1;
            Time_var(jj+100) = exp_no_treatment(ii,1);
        end
        tmp = tmp + k;
    end
end

simulation_survival_rate = ones(1, 200) * 100;
dead_num = zeros(200);
threshold = 0.695;

for t = 2:1:60
    k = 0;
    for person=1:1:100
        if simulation_no_treatment.C(person,t+1) > threshold
            k = k + 1;
            Event_var(person) = 1;
            Time_var(person) = t;
            simulation_no_treatment.C(person,:)=zeros(1,201);
        end
    end
    dead_num(t)=dead_num(t-1) + k;
    simulation_survival_rate(t) = 100 - dead_num(t);
end
%fig
Time_var = Time_var';
Event_var = Event_var';
[p2,~,~] = MatSurv(Time_var, Event_var, Group_var);


figure;
hold on
set(gca,'FontSize',15);
plot(1:200,simulation_survival_rate,'color','#34B34A',LineWidth=2.5);hold on;
plot(exp_no_treatment(:,1),exp_no_treatment(:,2),'color','#3855A5',LineWidth=2.5);hold on;
set(gca,'FontSize',14,'FontWeight','bold','fontname','Arial');
xlabel('Time (Days)','FontWeight','Bold','FontSize',16,'fontname','Arial');
ylabel('Percent Survival(%)','FontWeight','Bold','FontSize',16,'fontname','Arial');
title('Prediction','FontWeight','Bold','FontSize',16,'fontname','Arial');
% legend('No treatment','CSF1R-I','fontname','Arial','box','off');
ylim([0 100])
xlim([0 60])
text(0.63, 0.75, ['p-value = ', num2str(p2 ,'%.3g')], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',12,'FontWeight','bold','fontname','Arial');
box on;
hold off;

threshold = 0.695;
dead_num = zeros(2,4);
dead_num(1,1:4)=[0 66 34 0]; % exp
simulation_no_treatment=load('stochastic_No.mat');
for ii = 1:1:3
    k = 0;
    for person = 1:100
        for jj = 1:1:20
            if(simulation_no_treatment.C(person,(ii-1)*20+jj) > threshold)
                simulation_no_treatment.C(person,:)=zeros(1,201);
                k = k+1;
            end
        end
    end
    dead_num(2,ii+1)= k;  
end

figure;
set(gca,'FontSize',15); 
scatter(0:1:3,dead_num(2,1:4),'red','filled');
plot(0:1:3,dead_num(2,1:4),'Color','red','LineWidth',2);
xticks([0 1 2 3])
xticklabels({'0','20','40','60'}) 
xlabel('Time (Days)','FontWeight','Bold','FontSize',20);
xlim([0 3]);
ylim([0 100]);
hold on;
scatter(0:1:3,dead_num(1,1:4),'blue','filled');
plot(0:1:3,dead_num(1,1:4),'Color','blue','LineWidth',2);
legend('Prediction','','Experiment','');

box off;
hold off;