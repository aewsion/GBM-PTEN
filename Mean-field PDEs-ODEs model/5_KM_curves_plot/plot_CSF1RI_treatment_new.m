clc
close all
addpath('MatSurv-master');

simulation_CSF1RI_treatment=load('stochastic_CSF1R-I.mat');
exp_CSF1RI_treatment=load('Drug_Dataset.csv');

Time_var=140*ones(1,190);
Event_var=zeros(1,190);
Group_var = repmat({'prediction'}, 190, 1); % 前100 prediction,后90 experience。
%Group_var
for i=101:1:190
    Group_var(i) = {'experience'};
end

ll = length(exp_CSF1RI_treatment);
tmp = 1;
for ii=1:1:(ll-1)
    k = exp_CSF1RI_treatment(ii,3) - exp_CSF1RI_treatment(ii+1,3);
    if(k > 0)
        for jj = tmp:1:(tmp+k-1)
            Event_var(jj+100) = 1;
            Time_var(jj+100) = exp_CSF1RI_treatment(ii,1);
        end
        tmp = tmp + k;
    end
end

simulation_survival_rate = ones(1, 200) * 100;
dead_num = zeros(200);
threshold = 0.695;

for t = 2:1:140
    k = 0;
    for person=1:1:100
        if simulation_CSF1RI_treatment.C(person,t+1) > threshold
            k = k + 1;
            Event_var(person) = 1;
            Time_var(person) = t;
            simulation_CSF1RI_treatment.C(person,:)=zeros(1,201);
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
plot(1:140,simulation_survival_rate(1:140),'color','#34B34A',LineWidth=2.5);hold on;
plot(exp_CSF1RI_treatment(:,1),exp_CSF1RI_treatment(:,2),'color','#3855A5',LineWidth=2.5);hold on;
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

dead_n = zeros(1,9600);