addpath('.\MatSurv-master');
clear("Group_var");
Group_var = repmat({'Prediction'}, 200, 1);%repmat重复副本
for ii=101:1:200
    Group_var(ii) = {'Basal'};
end
%prediction 100; Experiment 90 (science 2016)
clear("Event_var");
Event_var = zeros(1,200);

clear("Time_var");
Time_var = zeros(1,200)+200;%200*ones(1,200)

Experiment_num = load(".\Drug_Dataset.csv");%时间，存活率，存活人数


Experiment_num = Experiment_num';%3*56
Experiment_num(1,:) = Experiment_num(1,:) * 48;
ll = length(Experiment_num);%56

clear("data_cc");
%data_cc = load("D:\plottest\Survival\2023_05_10\test24\k32=14.3\k4=5.7\cc_num.csv");
%data_cc = load("C:\Users\85456\Desktop\basal_figure\cc_num.csv");
data_cc = load(".\cc_num.csv");
data_cc = data_cc(:,0*48+1:200*48)/2154;
data_c = data_cc;
data_m = data_cc;
data_basal = load("D:\plottest\Survival\2023_05_10\test25\k32=14.3\k4=5.7\cc_num.csv");
data_n = data_basal(:,0*48+1:200*48)/2154;
clear("data_p");
clear("Survival_Rate");
clear("Basal_Rate");
S = sum(data_cc(:,48*10+1:48*50));%100个病人10天到50天的总和
[m,z] = min(S);
[M,Z] = max(S);
z = z + 48*10;
Z = Z + 48*10;
threshold = zeros(1,100);
for i = 1:100
    %R(i) = unifrnd(-0.05,0.05);
    %threshold(i) = 0.95+R(i);
    threshold(i) = 0.896;
end
c=jet(100);%颜色
figure;
pp = 0;
for i=1:100
    if(data_cc(i,140*48)==0)
        plot(data_cc(i,1:48*200),'Color','black','LineWidth',0.0001);
        
    else
        pp = pp + 1;
        plot(data_cc(i,1:48*200),'Color',c(i,:),'LineWidth',0.01);
        data_p(pp,:) = data_cc(i,1:48*200);
    end
    xticks([0 20*48 40*48 60*48 80*48 100*48 120*48 140*48 160*48 180*48 200*48])
    xticklabels({'0','20','40','60','80','100','120','140','160','180','200'}) 
    yticks([0  0.2  0.4  0.6  0.8  1.0])
    yticklabels({'0','0.2','0.4','0.6','0.8','1.0'})
    %yticks([0 215 430 646 861 1077 1292 1508 1723])
    %yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'})
    %yticks([45 90 135 180 225 270 315 360 405])
    %yticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'})    
    %xticks([0 2400 4800 7200 9600 12000])
    %xticklabels({'0','50','100','150','200','250'})
    xlim([0 200*48]);
    hold on;
            
end

set(gca,'FontSize',15); 
title("All Samples" ,'Color','black','FontSize',20,'FontAngle','normal','FontWeight','Bold');
xlabel('Time (Days)','FontWeight','Bold','FontSize',20);
ylabel("Tumor Density" ,'Color','black','FontSize',20,'FontAngle','normal','FontWeight','Bold');
%title("Average concentration of EGF" ,'Color','black','FontSize',16,'FontAngle','normal')
box off;
hold off;

figure;
x = Z:1:48*200;%Z=100个病人10到50天density和最大的那个病人
y = sum(data_p(:,Z:48*200));
y = y/pp;
[p,s] = polyfit(x,y,5);
y1 = polyval(p,x);
[yfit,dy] = polyconf(p,x,s,'predopt','curve');
for i=1:pp
    if(mod(i,20) == 0)
        plot(data_p(i,:),'Color','#6495ED','LineWidth',0.01);
        %plot(data_cc(i,1:9600),'Color','#EEF183','LineWidth',1.2);%M1
        %plot(data_cc(i,1:9600),'Color','#D4B261','LineWidth',0.8);%M2
    xticks([0 20*48 40*48 60*48 80*48 100*48 120*48 140*48 160*48 180*48 200*48])
    xticklabels({'0','20','40','60','80','100','120','140','160','180','200'}) 
        yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])
        yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
        ylim([0.4 1.0])
        %plot(data_cc(i,1:9600),'Color','#6495ED','LineWidth',1);
        hold on;
        xlim([0 200*48]);
        
   end
    
end
hold on;
h1 = fill([x,fliplr(x)],[yfit-1*dy,fliplr(yfit+1*dy)],'red');
hold on;
set(h1,'edgealpha',0,'facealpha',1);
set(gca,'FontSize',15); 
title("Drug Resistance Samples" ,'Color','black','FontSize',20,'FontAngle','normal','FontWeight','Bold');
xlabel('Time (Days)','FontWeight','Bold','FontSize',20);
ylabel("Tumor Density" ,'Color','black','FontSize',20,'FontAngle','normal','FontWeight','Bold');
%title("Average concentration of EGF" ,'Color','black','FontSize',16,'FontAngle','normal');
box off;
hold off;


figure;
dead_num = zeros(4,11);
dead_num(4,1) = 1;



   
for ii = 1:1:10
    k = 0;
    kk = 1;
    for l = 1:ll
        if(Experiment_num(1,l)>= 48*20*ii && kk==1)
            dead_num(4,ii+1) = l;
            dead_num(3,ii+1) = Experiment_num(2,dead_num(4,ii))-Experiment_num(2,dead_num(4,ii+1));
            %disp(l);
            %disp(dead_num);
            kk = 0;
        end
    end
    for j = 1:100
   
        for jj = 1:1:48*20
            if((ii-1)*48*20+jj > z)
            if(data_c(j,(ii-1)*48*20+jj)>= threshold(j))
                data_c(j,:)=zeros(1,48*200);
                k = k+1;
            end
            end
        end
    end
    %dead_num(1,ii)= 25+50*(ii-1);
    dead_num(2,ii+1)= k;
    
    hold on;
end
%%
% for ii = 1:1:3
%     k = 0;
%     for j = 1:20
%         avg_cc = 0;
%         for jj = 1:1:48*20
%             %avg_cc = sum(data_c(j,1:(ii-1)*48*20+jj))/((ii-1)*48*20+1+jj);
%             if (ii-1)*48*20+jj>z
%                 avg_cc = sum(data_c(j,(ii-1)*48*20+jj-1*48:(ii-1)*48*20+jj))/(1*48);
%             end
%             
%             if(avg_cc>= threshold)
%                 data_c(j,:)=zeros(1,48*70);
%                 k = k+1;
%             end 
%         end
%     end
%     dead_num(1,ii)= 25+50*(ii-1);
%     dead_num(2,ii+1)= k*5;
%     
%     hold on;
% end
%%
set(gca,'FontSize',15); 
scatter(0:1:10,dead_num(2,1:11),'red','filled');
plot(0:1:10,dead_num(2,1:11),'Color','red','LineWidth',2);
xticks([0 1 2 3 4 5 6 7 8 9 10])
xticklabels({'0','20','40','60','80','100','120','140','160','180','200'}) 
xlabel('Time (Days)','FontWeight','Bold','FontSize',20);
ylabel("Death Frequency" ,'Color','black','FontSize',20,'FontAngle','normal','FontWeight','Bold');
xlim([0 10]);
ylim([0 100]);
dead_num(1,1:8)=[0 0 0 9 25 21 1 0];
scatter(0:1:10,dead_num(1,1:11),'blue','filled');
plot(0:1:10,dead_num(1,1:11),'Color','blue','LineWidth',2);
%scatter(0:1:7,dead_num(3,1:8),'black','filled');
%plot(0:1:7,dead_num(3,1:8),'Color','black','LineWidth',2);
legend('Prediction','','Experiment','');
box off;
hold off;



%Survival_Rate = zeros(48*140);
for t=0*48+1:z
    Survival_Rate(t)=100;
    Basal_Rate(t)=100;
end

dead_n = zeros(2,9600);
for t=z+1:200*48
    k1 = 0;
    k2 = 0;
    for i=1:100
        avg_cc = sum(data_m(i,t-1:t))/(2);
        %avg_cc
        if(data_m(i,t)>=threshold(j))
        %if(data_m(i,t)>=threshold && data_m(i,t-20*48)>=threshold)
            k1 = k1 + 1;
            data_m(i,:)=zeros(1,48*200);
            Event_var(i) = 1;
            Time_var(i) = t/48;
        end
        if(data_n(i,t)>=threshold(j))
            k2 = k2 + 1;
            data_n(i,:)=zeros(1,48*200);
            Event_var(i+100) = 1;
            Time_var(i+100) = t/48;
        end
    end
    dead_n(1,t)=dead_n(1,t-1)+k1;
    Survival_Rate(t)=100-dead_n(1,t);
    dead_n(2,t)=dead_n(2,t-1)+k2;
    Basal_Rate(t)=100-dead_n(2,t);
% if t>=2
% if Survival_Rate(t-1)==min(Survival_Rate)
%     Survival_Rate(t)=0;
% end 
% end
end
Time_var = Time_var';
Event_var = Event_var';
[p2,~,~] = MatSurv(Time_var, Event_var, Group_var);

figure;
hold on
set(gca,'FontSize',15); 
plot(1:48*200,Survival_Rate,'color','red',LineWidth=2);
plot(1:48*200,Basal_Rate,'color','black',LineWidth=2);
%plot(Experiment_num(1,:),Experiment_num(2,:),'Color','blue','LineWidth',2);
xlabel('Time (Days)','FontWeight','Bold','FontSize',20);
ylabel('Percent Survival','FontWeight','Bold','FontSize',20);
    
figure_FontSize=20;
xticks([0 20*48 40*48 60*48 80*48 100*48 120*48 140*48 160*48 180*48 200*48])
xticklabels({'0','20','40','60','80','100','120','140','160','180','200'}) 
xlim([0 200*48])
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
yticks([0 20 40 60 80 100]);
yticklabels({'0','20','40','60','80','100'});
ylim([0 100]);
%set(findobj('FontSize',14),'FontSize',figure_FontSize);
legend('Prediction','Basal');
text(0.75, 0.5, ['p-value = ', num2str(p2 ,'%.3g')], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
box off;
hold off;

