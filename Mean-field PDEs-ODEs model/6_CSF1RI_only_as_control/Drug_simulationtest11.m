 function pde_Paro_Ellip

clc
% close all
clear all

global C_T_mean  %整个坐标轴上肿瘤细胞的平均密度
global eta1    %药物的降解率
global D_d1    %药物的扩散率
global delt_t  %模拟时候的时间间隔
global L      %模拟时候横轴x的大小
global B_2  %这个是需要估计的值,是IL4分泌速率相关的
global A_0   %这个是需要估计的值,是耐药细胞因子初始值
global d_max         %药物边界值，等于3
global d_initial     %和药物初始的量相关的，初始量为d_initial*x，且需d_initial*L=d_max
global period        %这个是用于边界用药时候有一个周期，例如7天一个周期
global inter         %阶梯函数的换化间隔，例如用一天把药物量为0渐变成药物量为6，一天就是这个inter
global period1       %真实的的时间周期，由period×delt_t得到
global inter1        %真实的的时间间隔，由inter×delt_t得到
global TimeLength    %总的时间
global m             %边界条件三设计的药物的量
global Ptest1
global EGFR
global ERK
global IGF1R
global AKT
global EGF
global IGF1

load('D:\论文\老板给的文章\老板给的2013年文章自己改\Model 3 - 改进模型最终2\canshu111.mat','p');
Ptest1=p(10,1:22);
ExpData;
EGFR_0=pY1068_EGFR(1,1);
ERK_0=p_ERK(1,1);
IGF1R_0=p_IGF1Rbeta(1,1);
AKT_0=pS473_AKT(1,1);


m=3;%这个是第三种情况下，阶梯函数的最高值，代表喂药一段时间的量为m，一段时间的量为0
TimeLength=200;        %总时间
TimeLength1=200;       %画图时候坐标的大小
delt_t=1/91.5926;      %模拟时候的时间间隔
%delt_t是一天，因为3600*24/t_s= 0.0108=1/92.5926
period=20;           %小周期，例如20天用药量6，20天用药量0
period1=period*delt_t;  %真正周期
inter=1;             %间隔
inter1=inter*delt_t;   %真正间隔
x_gap=0.1;             %模拟时候坐标轴的间隔
L=100;         %模拟时候坐标轴的总长度
x = 0:x_gap:L;           %肿瘤生长的空间，一维的
%x1 = 5:1:105
% x=0:0.01:1;

Parameters_nondimensional_3_gai_jin_fangcheng_3
A_0=A_0;
t = 0:delt_t:delt_t*TimeLength;
eta1=eta_d;
D_d1=D_d;%药物的扩散率
d_max=3;%药物边界的量
d_initial=d_max/L;%和药物初始的量相关的，初始量为d_initial*x，且需d_initial*L=d_max
B_2=B_2;
 
xhinter1=20;
xhinter2=100;
yhinter1=20;          %画图间隔
yhinter2=200;         %画图间隔
Conditions;
ExpData;

     EGF=Cond1(ii,1);
     IGF1=Cond1(ii,2);
     EGFR_I=Cond1(ii,3);
     IGF1R_I=Cond1(ii,4);
  


options=odeset('reltol',1e-3);%精度函数不能少
m0 = 2;%用matlab解需要设置的参数
sol = pdepe(m0,@pdex4pde,@pdex4ic,@pdex4bc,x,t,options);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);
u5 = sol(:,:,5);
u6 = sol(:,:,6);
U=[u1;u2;u3;u4;u5;u6];
 C_T_mean=sum(u1(:,:),2)/1001; 
 
 
%    figure,
%    hold on
%    plot(0:TimeLength,Drug_mean,'m-','LineWidth',3);
%    label=num2str(0:yhinter1:TimeLength);
%    set(gca,'xtick',[0:yhinter1:TimeLength])
%    set(gca,'xticklabel',0:yhinter1:TimeLength1)
%    xlabel('Days','FontWeight','Bold','FontSize',18);ylabel('Drug_mean','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
              
% figure
% %%%%% surf(x,t,u1)
% pcolor(x,t,A1)
% title('Drugt1')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
%    figure
% % % % surf(x,t,u1)
% pcolor(x,t,A2)
% title('Drugt2')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);

% figure
% % % % surf(x,t,u1)
% pcolor(x,t,A3)
% title('Drugt3')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % surf(x,t,u1)
% pcolor(x,t,A4)
% title('Drugt4')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);

% figure
% % % surf(x,t,u1)
% pcolor(x,t,A5)
% title('Drug-3')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% % 
% figure
% % % surf(x,t,u1)
% pcolor(x,t,A6)
% title('IL4_3')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % surf(x,t,u1)
% pcolor(x,t,A7)
% title('IL4_1')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);
% 
% figure
% % % surf(x,t,u1)
% pcolor(x,t,A8)
% title('IL4_2')
% colorbar; 
% shading interp;
%    set(gca,'xtick',[0:20:100]); 
%    set(gca,'xticklabel',0:20:100);
%    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
%    set(gca,'yticklabel',0:yhinter1:TimeLength1)
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
%    set(gca,'FontWeight','Bold','FontSize',12);


  
   
%     axis([0 200 0 1])

figure
% surf(x,t,u1)
pcolor(x,t,u1)
title('Tumor')
colorbar; shading interp;
   set(gca,'xtick',[0:20:100]); 
   set(gca,'xticklabel',0:20:100);
   set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); 
   set(gca,'yticklabel',0:yhinter1:TimeLength1)
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',12);


figure
% surf(x,t,u2)
pcolor(x,t,u2)
title('Macophagy_1')
colorbar; shading interp;
   set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); set(gca,'yticklabel',0:yhinter1:TimeLength1)
   set(gca,'xtick',[0:20:100]); set(gca,'xticklabel',0:20:100);
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',12);


figure
% surf(x,t,u3)
pcolor(x,t,u3)
title('Macophagy_2')
colorbar; shading interp;
   set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); set(gca,'yticklabel',0:yhinter1:TimeLength1)
   set(gca,'xtick',[0:20:100]); set(gca,'xticklabel',0:20:100);
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',12);


figure
% surf(x,t,u4)
pcolor(x,t,u4)
title('CSF1')
colorbar; shading interp;
   set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); set(gca,'yticklabel',0:yhinter1:TimeLength1)
   set(gca,'xtick',[0:20:100]); set(gca,'xticklabel',0:20:100);
   xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',12);


 figure
 % surf(x,t,u5)
 pcolor(x,t,u5)
 title('EGF')
 colorbar; shading interp;
    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); set(gca,'yticklabel',0:yhinter1:TimeLength1)
    set(gca,'xtick',[0:20:100]); set(gca,'xticklabel',0:20:100);
    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
    set(gca,'FontWeight','Bold','FontSize',12);
    
 figure
 % surf(x,t,u5)
 pcolor(x,t,u6)
 title('IGF1')
 colorbar; shading interp;
    set(gca,'ytick',0:delt_t*yhinter1:delt_t*TimeLength); set(gca,'yticklabel',0:yhinter1:TimeLength1)
    set(gca,'xtick',[0:20:100]); set(gca,'xticklabel',0:20:100);
    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
    set(gca,'FontWeight','Bold','FontSize',12);
end
% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx)
global d_max
global B_2
global Ptest1
global EGFR
global ERK
global IGF1R
global AKT
ExpData;
q=0;
qq=0;
qqq=0;
% global drug
Parameters_nondimensional_3_gai_jin_fangcheng_3
c = [1; 1; 1; 1; 1; 1];  
% Conditions_ODES;

 %size(D 0,0,0,1/e,0;uDx)
f = [0,0,0,0;
    0,0,0,0;
    0,0,0,0;
    0,0,0,0]*DuDx; 
t
%调试
% %      R1=r_T*(1+a_r*u(5)./(K1+u(5))+a_r*u(6)./(K4+u(6)));  % growth rate
%      R1=r_T*(1+a_r*ERK./(K5+ERK)+a_r*AKT./(K6+AKT));
%      D1=d_T*(1+k_d*u(2));  % death rate    
% %      R3=k_M*u(3)-u(2).*(a_CSF1*u(4)./(K21+K22*drug+u(4))+a_A*IL4/(K3+IL4));
%      R3=k_M*u(3)-u(2).*(a_CSF1*u(4)./(K21+K22*drug+u(4))+a_A*IL4/(K3+IL4));
%      bs=u(3)*IL4-u(6);
%     if bs<0
%         bs=0;
%     end
    R1=V1*EGF/(K11+EGF)/(1+y(2)/K12)/(1+EGFR_I/K13)*(1-y(1))-d1*y(1);                    %EGFR激活的
    R2=a*(1+V21*y(1)^n/(K21^n+y(1)^n))*(1+V22*y(3)/(K22+y(3)))/(1+y(4)/K23)*(1-y(2))-d2*y(2);    %ERK激活的
    R3=V3*IGF1/(K31+IGF1)/(1+y(2)/K32)/(1+IGF1R_I/K33)*(1-y(3))-d3*y(3);                %IGFR激活的
    R4=V4*y(1)/(K41+y(1))*y(3)/(K42+y(3))*(1-y(4))-d4*y(4);                                 %AKT激活的

% 新的
 s = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)-u(4));1/epsilon*B_3*(u(3).*u(4)./(K21+K22*drug+u(4))-u(5));1/epsilon*bs]; 
end
% --------------------------------------------------------------
function u0 = pdex4ic(x)
global L
d_max=3;
d_initial=d_max/L;
for iii=1:length(x)
    if iii<=11
        C0(iii)=0.5;
    else
        C0(iii)=0;
    end
end
% u0 = [0.5.*Heaviside(11-x).*(1+0.1.*randn(1,length(x))); 100/100.*(1+0.1.*randn(1,length(x)))+0.0; 1.*(1+0.01.*randn(1,length(x))); 1.*(1+0.01.*randn(1,length(x))); d_max*x/100];    
u0 = [0.5.*(1+0.1.*rand(1,length(x))); x/100.*(1+0.2.*randn(1,length(x)))+2; x/100.*(1+0.2.*randn(1,length(x)))+2;0.2.*(1+0*0.01.*rand(1,length(x)));0.85.*(1+0*0.01.*rand(1,length(x)));0]; 
% u0 = [0.5.*(1+0.1.*rand(1,length(x))); x/100.*(1+0.1.*randn(1,length(x)))+1; x/100.*(1+0.1.*randn(1,length(x)))+1;1.*(1+0*0.01.*rand(1,length(x))); 1.*(1+0*0.01.*rand(1,length(x)))]; 
%u0 = [0.5.*(1+0.1.*rand(1,length(x))); x/100.*(1+0.1.*randn(1,length(x)))+1; 1.*(1+0*0.01.*rand(1,length(x))); 1.*(1+0*0.01.*rand(1,length(x))); d_max*x/100]; 
%u0 = [0.5.*(1+0.1.*rand(1,length(x))); x/100.*(1+0.1.*randn(1,length(x)))+1; 1.*(1+0*0.01.*rand(1,length(x))); 1.*(1+0*0.01.*rand(1,length(x))); d_max*x/100]; 
%u0 = [0.5.*(1+0.1.*rand(1,length(x))); x/100.*(1+0.1.*randn(1,length(x)))+1; 1.*(1+0*0.01.*rand(1,length(x))); 1.*(1+0*0.01.*rand(1,length(x))); d_max*x/100]; 
end


% --------------------------------------------------------------

function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
global d_max
global delt_t
global period1
global inter1
pl = [0; 0; 0; 0; 0; 0];                               
ql = [1; 1; 1; 1; 1; 1];
pr = [0; 0; 0; 0; 0; 0];
%pr = [0; 0; 0; 0; ur(5)-dnt];
%pr = [0; 0; 0; 0; ur(5)-d_max/2*(sin((t/delt_t)*pi/100)+2)];
%pr = [0; 0; 0; 0; ur(5)-(3-rem(floor(t/(delt_t*30)),3))];%忘记了是哪种情况了
qr = [1; 1; 1; 1; 1; 1];
end

% --------------------------------------------------------------
%%以下函数应该是边界条件分别为d11=d_max=3和d22=d_max/2*(sin(w*t)+2)所计算出来药物在x坐标轴上的浓度
function [drugt1,drugt2] = picture(x,t)
global d_max
global eta1
global D_d1
global L   %x的最大值
global delt_t
global TimeLength
global d_initial
w=pi/(TimeLength*delt_t/20);
   d11=d_max;
   d110=d_max;
   d22=d_max/2*(sin(w*t)+2);
drugt1=0;
drugt2=0;
drugt11=0;
drugt21=0;
for n=1:1:100
    Tn1=((2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+  d110*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t)))*sin(n*pi*x/L);
    drugt11=drugt11+Tn1; 
    Tn2=(((2*L^2*(-1)^(n+1)/(n*pi)+4*L^2*(-1)^n/(n*pi)^3-4*L^2/(n*pi)^3)*d_initial+2*L*(-1)^n/(n*pi)*d_max)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t))/(n*pi*D_d1*(n*pi/L)^2)+d_max*w*2*L*(-1)^n*(w*sin(w*t)+D_d1*(n*pi/L)^2*(cos(w*t)-exp(-D_d1*(n*pi/L)^2*t)))/(2*n*pi*(w^2+D_d1^2*(n*pi/L)^4)))*sin(n*pi*x/L);
    drugt21=drugt21+Tn2;
end
drugt1=1/x*drugt11+d11;
drugt2=1/x*drugt21+d22;
end

%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用的sin函数来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分
function q = Tn22(x,t,n)
global D_d1;
global eta1;
global delt_t
global L
global period1
global inter1
q=0;
q1=0;
q2=0;
d_0=6;          %药物的边界值
d_initial=d_0/L;%药物的初始条件为d_initial*r，且需要d_initial*L=d_0
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
for i=0:1:(d2-1)
    q1=6*L/(n*pi)*(-1)^n*(-exp(D_d1*(n*pi/L)^2*(i*2+1)*period1-D_d1*(n*pi/L)^2*t)-exp(D_d1*(n*pi/L)^2*((i*2+1)*period1-inter1)-D_d1*(n*pi/L)^2*t))*L^4/(L^4+inter1^2*D_d1^2*n^4*pi^2);
    q2=q2+q1;  
    q1=6*L/(n*pi)*(-1)^n*(exp(D_d1*(n*pi/L)^2*(i*2+2)*period1-D_d1*(n*pi/L)^2*t)+exp(D_d1*(n*pi/L)^2*((i*2+2)*period1-inter1)-D_d1*(n*pi/L)^2*t))*L^4/(L^4+inter1^2*D_d1^2*n^4*pi^2);
    q2=q2+q1;
end
if d3==0 && t0>=(period1-inter1)
    q1=6*L/(n*pi)*(-1)^n*(sin((t+3*inter1/2-(d1+1)*period1)*pi/inter1)*exp(D_d1*(n*pi/L)^2*t-D_d1*(n*pi/L)^2*t)-exp(D_d1*(n*pi/L)^2*((d1+1)*period1-inter1)-D_d1*(n*pi/L)^2*t)+inter1/pi*D_d1*(n*pi/L)^2*cos((t+3*inter1/2-(d1+1)*period1)*pi/inter1)*exp(D_d1*(n*pi/L)^2*t-D_d1*(n*pi/L)^2*t))*L^4/(L^4+inter1^2*D_d1^2*n^4*pi^2);
    q2=q2+q1; 
end
if d3==1 && t0<(period1-inter1)
    q1=6*L/(n*pi)*(-1)^n*(-exp(D_d1*(n*pi/L)^2*(d2*2+1)*period1-D_d1*(n*pi/L)^2*t)-exp(D_d1*(n*pi/L)^2*((d2*2+1)*period1-inter1)-D_d1*(n*pi/L)^2*t))*L^4/(L^4+inter1^2*D_d1^2*n^4*pi^2);
    q2=q2+q1; 
end 
if d3==1 && t0>=(period1-inter1)
    q1=6*L/(n*pi)*(-1)^n*(-exp(D_d1*(n*pi/L)^2*(d2*2+1)*period1-D_d1*(n*pi/L)^2*t)-exp(D_d1*(n*pi/L)^2*((d2*2+1)*period1-inter1)-D_d1*(n*pi/L)^2*t))*L^4/(L^4+inter1^2*D_d1^2*n^4*pi^2);
    q2=q2+q1; 
    q1=6*L/(n*pi)*(-1)^n*(sin((t+inter1/2-(d1+1)*period1)*pi/inter1)*exp(D_d1*(n*pi/L)^2*t-D_d1*(n*pi/L)^2*t)+exp(D_d1*(n*pi/L)^2*((d1+1)*period1-inter1)-D_d1*(n*pi/L)^2*t)+inter1/pi*D_d1*(n*pi/L)^2*cos((t+inter1/2-(d1+1)*period1)*pi/(inter1))*exp(D_d1*(n*pi/L)^2*t-D_d1*(n*pi/L)^2*t))*L^4/(L^4+inter1^2*D_d1^2*n^4*pi^2);
    q2=q2+q1; 
end 
q=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t))+q2;
end
%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用的sin函数来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q1 = Drug_22(x,t)  
global L
global period1
global inter1
q=0;
dnt=0;
%al=floor(t/(delt_t*period));

d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
if d3==0 && t0<=(period1-inter1)
    dnt=6;
end
if d3==0 && t0>(period1-inter1)
    dnt=3+3*sin((t+3*inter1/2-(d1+1)*period1)*pi/inter1);
    
end
if d3==1 && t0<=(period1-inter1)
    dnt=0;
end
if d3==1 && t0>(period1-inter1)
    dnt=3+3*sin((t+inter1/2-(d1+1)*period1)*pi/inter1);
    
end
        
    for n=1:1:20
    q=q+Tn22(x,t,n)*sin(n*pi*x/L);  
    end

q1=q/x+dnt;
end

%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分
function q2=air(x1,x2,y1,y2,t1,t2,t,n)
global D_d1
global L
x=0;
q1=0;
q2=0;
x=t1;
d11x=2*y1*((x-x2)^2+(x2-x1)*(x-x2)+2*(x-x1)*(x-x2))/(x2-x1)^3+2*y2*((x-x1)^2+(x1-x2)*(x-x1)+2*(x-x2)*(x-x1))/(x1-x2)^3;
d12x=2*y1*(4*(x-x2)+(x2-x1)+2*(x-x1))/(x2-x1)^3+2*y2*(4*(x-x1)+(x1-x2)+2*(x-x2))/(x1-x2)^3;
d13x=12*y1/(x2-x1)^3+12*y2/(x1-x2)^3;
x=t2;
d21x=2*y1*((x-x2)^2+(x2-x1)*(x-x2)+2*(x-x1)*(x-x2))/(x2-x1)^3+2*y2*((x-x1)^2+(x1-x2)*(x-x1)+2*(x-x2)*(x-x1))/(x1-x2)^3;
d22x=2*y1*(4*(x-x2)+(x2-x1)+2*(x-x1))/(x2-x1)^3+2*y2*(4*(x-x1)+(x1-x2)+2*(x-x2))/(x1-x2)^3;
d23x=12*y1/(x2-x1)^3+12*y2/(x1-x2)^3;
q1=(d21x*exp(D_d1*(n*pi/L)^2*(t2-t))-d11x*exp(D_d1*(n*pi/L)^2*(t1-t)))*2*L/(n*pi)*(-1)^n/(D_d1*(n*pi/L)^2);
q2=q2+q1;
q1=(d22x*exp(D_d1*(n*pi/L)^2*(t2-t))-d12x*exp(D_d1*(n*pi/L)^2*(t1-t)))*2*L/(n*pi)*(-1)^n/(D_d1^2*(n*pi/L)^4);
q2=q2-q1;
q1=(d23x*exp(D_d1*(n*pi/L)^2*(t2-t))-d13x*exp(D_d1*(n*pi/L)^2*(t1-t)))*2*L/(n*pi)*(-1)^n/(D_d1^3*(n*pi/L)^6);
q2=q2+q1;
end
%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q = Tn23(t,n)
global D_d1;
global eta1;
global L
global period1
global inter1
q=0;
q1=0;
q2=0;
d_0=6;
%初始条件为d_initial*r，应该是d_0除以L
d_initial=d_0/L;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
for i=0:1:(d2-1)
    x1=(2*i+1)*period1-inter1;
    x2=(2*i+1)*period1;
    y1=6;
    y2=0;
    t1=x1;
    t2=x2;
    q1=air(x1,x2,y1,y2,t1,t2,t,n);
    q2=q2+q1;
    x1=(2*i+2)*period1-inter1;
    x2=(2*i+2)*period1;
    y1=0;
    y2=6;
    t1=x1;
    t2=x2;
    q1=air(x1,x2,y1,y2,t1,t2,t,n);
    q2=q2+q1;
end
if d3==0 && t0>=(period1-inter1)
    x1=(d1+1)*period1-inter1;
    x2=(d1+1)*period1;
    y1=6;
    y2=0;
    t1=x1;
    t2=t;
    q1=air(x1,x2,y1,y2,t1,t2,t,n);
    q2=q2+q1; 
end
if d3==1 && t0<(period1-inter1)
    x1=d1*period1-inter1;
    x2=d1*period1;
    y1=6;
    y2=0;
    t1=x1;
    t2=x2;
    q1=air(x1,x2,y1,y2,t1,t2,t,n);
    q2=q2+q1;
end 
if d3==1 && t0>=(period1-inter1)
    x1=d1*period1-inter1;
    x2=d1*period1;
    y1=6;
    y2=0;
    t1=x1;
    t2=x2;
    q1=air(x1,x2,y1,y2,t1,t2,t,n);
    q2=q2+q1;
    x1=(d1+1)*period1-inter1;
    x2=(d1+1)*period1;
    y1=0;
    y2=6;
    t1=x1;
    t2=t;
    q1=air(x1,x2,y1,y2,t1,t2,t,n);
    q2=q2+q1;
end 
q=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t))+q2;
end
%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q1 = Drug_23(x,t) 
global L
global period1
global inter1
q=0;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
xk=(d1+1)*period1-inter1;
xk1=(d1+1)*period1;
if d3==0 && t0<=(period1-inter1)
    dnt=6;
end
if d3==0 && t0>(period1-inter1)
    yk=6;
    %yk1=0;
    dnt=(1+2*(t-xk)/(xk1-xk))*((t-xk1)/(xk-xk1))^2*yk+0+0+0;   
end
if d3==1 && t0<=(period1-inter1)
    dnt=0;
end
if d3==1 && t0>(period1-inter1)
    %yk=0;
    yk1=6;
    dnt=0+(1+2*(t-xk1)/(xk-xk1))*((t-xk)/(xk1-xk))^2*yk1+0+0;
    
end
    for n=1:1:20
    q=q+Tn23(t,n)*sin(n*pi*x/L);  
    end
q1=q/x+dnt;
end

  
%下面是计算在位置x和时间t时由药物CSF1R_I产生的药物产生的IL4的浓度

%%药物边界条件为常量3的时候产生的IL4的量
%第一部分
function q = TnIL4_1(x,t,n)

global D_d1;
global eta1;
global L
q=0;
d_0=3;%药物右边界条件，持续输入
d_initial=d_0/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
% q2=0;%因为又边界输入为常数3，导数值为0；
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/(D_d1*(n*pi/L)^2)+eta1*2*(-1)^n*t/(D_d1*(n*pi/L)^3)+eta1*2*(-1)^n*(exp(-D_d1*(n*pi/L)^2*t)-1)/(D_d1^2*(n*pi/L)^5);

end
%第二部分
function q = IL4_1(x,t)
global A_0
global L
global B_2  %IL4的分泌速率 
q0=0;
d_0=3;%初始右边界药物的条件

for n=1:1:100
    q0=q0+sin(n*pi/L*x)*TnIL4_1(x,t,n);
end
q=A_0+B_2/x*q0+B_2*d_0*t;
end

%%药物边界条件为d22=d_max/2*(sin(w*t)+2)=3/2*(sin(w*t)+2)的时候产生的IL4的量
%第一部分
function q5 = TnIL4_2(t,n)

global D_d1;
global eta1;
global L
global delt_t
global TimeLength
w=pi/(TimeLength*delt_t/20);
q=0;
d_0=3;%药物右边界条件，持续输入
d_initial=d_0/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
% q2=0;%因为又边界输入为常数3，导数值为0；
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q1=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/(D_d1*(n*pi/L)^2);
q2=eta1*2*(-1)^n*t/(D_d1*(n*pi/L)^3);
q3=eta1*2*(-1)^n/(D_d1^2*(n*pi/L)^5)*(exp(-D_d1*(n*pi/L)^2*t)-1);
q4=d_0*w*L/n/pi*(-1)^n*(D_d1*(n*pi/L)^2*sin(w*t)/w+exp(-D_d1*(n*pi/L)^2*t)-cos(w*t))/(w^2+D_d1^2*(n*pi/L)^4);
q5=q1+q2+q3+q4;
end
%第二部分
function q = IL4_2(x,t)
global A_0
global L
global B_2  %IL4的分泌速率 
global TimeLength
global delt_t
q0=0;
d_0=3;%初始右边界药物的条件
w=pi/(TimeLength*delt_t/20);
for n=1:1:100
    q0=q0+sin(n*pi*x/L)*TnIL4_2(t,n);
end
q=A_0+B_2/x*q0+B_2*d_0/2/w*(1-cos(w*t))+B_2*d_0*t;
end

%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分
%积分函数一，对应h1(t,t1)  t1是to
function q4=h_1_ceshiair1(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*exp(-D_d1*(n*pi/L)^2*(t-t1))*(2*m/d^3*(t1-(k+1/2)*T)^2+2*m/d^2*(1+2/d*(t1-(k+1/2)*T+d))*(t1-(k+1/2)*T));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*exp(-D_d1*(n*pi/L)^2*(t-t1))*(8*m/d^3*(t1-(k+1/2)*T)+2*m/d^2*(1+2/d*(t1-(k+1/2)*T+d)));
q3=24*m*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(t-t1));
q4=q1+q2+q3;
end
%积分函数二，对应h1(t,t1)
function q4=h_2_ceshiair2(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*exp(-D_d1*(n*pi/L)^2*(t-t1))*(-2*m/d^3*(t1-(k+1)*T+d)^2+2*m/d^2*(1-2/d*(t1-(k+1)*T))*(t1-(k+1)*T+d));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*exp(-D_d1*(n*pi/L)^2*(t-t1))*(-8*m/d^3*(t1-(k+1)*T+d)+2*m/d^2*(1-2/d*(t1-(k+1)*T)));
q3=-24*m*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(t-t1));
q4=q1+q2+q3;
end
%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q = Tn34(t,n)
global D_d1;
global eta1;
global L
global period1
global inter1
global m
T=2*period1;
d=inter1;
q=0;
q1=0;
q2=0;
q3=0;
q4=0;
d_0=m;%初始条件为d_initial*r，应该要除以L
d_initial=d_0/L;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
for i=0:1:(d2-1)
    t1=(i+1/2)*T-d;
    q1=h_1_ceshiair1(i,t1,t,n);
    t1=(i+1/2)*T;
    q2=h_1_ceshiair1(i,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(i+1)*T-d;
    q1=h_2_ceshiair2(i,t1,t,n);
    t1=(i+1)*T;
    q2=h_2_ceshiair2(i,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end
if d3==0 && t0>=(period1-inter1)
    t1=(d2+1/2)*T-d;
    q1=h_1_ceshiair1(d2,t1,t,n);
    t1=t;
    q2=h_1_ceshiair1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==1 && t0<(period1-inter1)
    t1=(d2+1/2)*T-d;
    q1=h_1_ceshiair1(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_1_ceshiair1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
end 
if d3==1 && t0>=(period1-inter1)
    t1=(d2+1/2)*T-d;
    q1=h_1_ceshiair1(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_1_ceshiair1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1)*T-d;
    q1=h_2_ceshiair2(d2,t1,t,n);
    t1=t;
    q2=h_2_ceshiair2(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end 
q=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t))+q4;
end
%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q = Drug_34(x,t) 
global L
global period1
global inter1
global m
q0=0;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
xk=(d1+1)*period1-inter1;
xk1=(d1+1)*period1;
if d3==0 && t0<=(period1-inter1)
    dnt=m;
end
if d3==0 && t0>(period1-inter1)
    yk=m;
    %yk1=0;
    dnt=(1+2*(t-xk)/(xk1-xk))*((t-xk1)/(xk-xk1))^2*yk+0+0+0;   
end
if d3==1 && t0<=(period1-inter1)
    dnt=0;
end
if d3==1 && t0>(period1-inter1)
    %yk=0;
    yk1=m;
    dnt=0+(1+2*(t-xk1)/(xk-xk1))*((t-xk)/(xk1-xk))^2*yk1+0+0;
    
end
    for n=1:1:20
    q0=q0+Tn34(t,n)*sin(n*pi*x/L);  
    end
q=q0/x+dnt;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数三
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%第一部分
function q6=h_3_ceshiair3(k,s,n)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*(2*m/3/d^3*(s-(k+1/2)*T)^3+2*m/d^2*(s^2/2+2/d*(s^3/3-((k+1/2)*T-d)/2*s^2))-2*m*(k+1/2)*T/d^2*(s+1/d*(s-(k+1/2)*T+d)^2));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*(4*m/d^3*(s-(k+1/2)*T)^2+2*m/d^2*(s+1/d*(s-(k+1/2)*T+d)^2));
q3=24*m*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*s;
q4=12*m*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q5=24*m*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q6=q1+q2+q3+q4+q5;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数四
function q6=h_4_ceshiair4(k,s,n)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*(-2*m/3/d^3*(s-(k+1)*T+d)^3+2*m/d^2*(s^2/2-2/d*(s^3/3-(k+1)*T/2*s^2))-2*m*((k+1)*T-d)/d^2*(s-1/d*(s-(k+1)*T)^2));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*(-4*m/d^3*(s-(k+1)*T+d)^2+2*m/d^2*(s-1/d*(s-(k+1)*T)^2));
q3=-24*m*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*s;
q4=-12*m*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q5=-24*m*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q6=q1+q2+q3+q4+q5;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数五
%对应的是对gk1函数的积分
%积分函数五是因为fn包含边界函数d'(t)是阶梯函数，分段积分完成了后是一个定的关于t的积分，再对t积分后和积分函数三不一样，然后才是IL的值的一部分
function q5=h_5_ceshiair5(k,s,n)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=12*m*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T));
q2=12*m*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q3=-24*m*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T));
q4=24*m*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数六
%对应的是对gk2函数的积分
function q5=h_6_ceshiair6(k,s,n)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=-12*m*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T));
q2=-12*m*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q3=24*m*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T));
q4=-24*m*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数七
%对应的是对d(t)阶梯函数第一部分的积分
function q1=h_7_ceshiair7(s)
global m
q1=m*s;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数八
%对应的是对d(t)阶梯函数第二部分的积分
function q5=h_8_ceshiair8(k,s)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=m/d^2/2/d*s^4;
q2=m/d^2/3*(3-2/d*(k+1/2)*T-4/d*(k+1/2)*T)*s^3;
q3=m/d^2/2*(2/d*(k+1/2)^2*T^2-2*(3-2/d*(k+1/2)*T)*(k+1/2)*T)*s^2;
q4=m/d^2*(3-2/d*(k+1/2)*T)*(k+1/2)^2*T^2*s;
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算IL4，对药物的一部分进行积分，下面为积分函数八
%对应的是对d(t)阶梯函数第二部分的积分
function q5=h_9_ceshiair9(k,s)
global D_d1
global L
global inter1
global period1
global m
T=2*period1;
d=inter1;
q1=-m/d^2/2/d*s^4;
q2=m/d^2/3*(1+2/d*(k+1)*T-4/d*(d-(k+1)*T))*s^3;
q3=m/d^2/2*(-2/d*(d-(k+1)*T)^2+2*(1+2/d*(k+1)*T)*(d-(k+1)*T))*s^2;
q4=m/d^2*(1+2/d*(k+1)*T)*(d-(k+1)*T)^2*s;
q5=q1+q2+q3+q4;
end
%%药物边界条件为阶梯函数6和0的时候产生的IL4的量，中间用埃尔米插值连接
%第二部分
function q = TnIL4_3(t,n)

global D_d1;
global eta1;
global L
global period1  %period1为一个周期
global inter1
global m
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q=0;
q4=0;
d_0=m;%药物右边界条件，持续输入
d_initial=d_0/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%积分函数五是因为fn包含边界函数d'(t)是阶梯函数，分段积分完成了后是一个定的关于t的积分，再对t积分后和积分函数三不一样，然后才是IL的值的一部分
for i=0:1:(d2-1)
    s=(i+1/2)*T-d;
    q1=h_3_ceshiair3(i,s,n);%调用积分函数三
    s=(i+1/2)*T;
    q2=h_3_ceshiair3(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T;
    q1=h_5_ceshiair5(i,s,n);%调用积分函数五
    s=t;
    q2=h_5_ceshiair5(i,s,n);%调用积分函数五
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T-d;
    q1=h_4_ceshiair4(i,s,n);%调用积分函数四
    s=(i+1)*T;
    q2=h_4_ceshiair4(i,s,n);%调用积分函数四
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T;
    q1=h_6_ceshiair6(i,s,n);%调用积分函数六
    s=t;
    q2=h_6_ceshiair6(i,s,n);%调用积分函数六
    q3=q2-q1;
    q4=q4+q3;
end
% for i=0:1:(d2-1)
%     s=(i+1/2)*T-d;
%     q1=h_3_ceshiair3(i,s,n);%调用积分函数三
%     s=(i+1/2)*T;
%     q2=h_3_ceshiair3(i,s,n);%调用积分函数三
%     q3=q2-q1;
%     q4=q4+q3;
%     s=(i+1)*T-d;
%     q1=h_4_ceshiair4(i,s,n);%调用积分函数四
%     s=(i+1)*T;
%     q2=h_4_ceshiair4(i,s,n);%调用积分函数四
%     q3=q2-q1;
%     q4=q4+q3;
% end

%判断t最后走到第几个周期积分对应着对应h1(s),h2(s)的积分
if d3==0 && t0>=(period1-inter1)
    s=(d2+1/2)*T-d;
    q1=h_3_ceshiair3(d2,s,n);%调用积分函数三
    s=t;
    q2=h_3_ceshiair3(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==1 && t0<(period1-inter1)
    s=(d2+1/2)*T-d;
    q1=h_3_ceshiair3(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_3_ceshiair3(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T;
    q1=h_5_ceshiair5(d2,s,n);%调用积分函数五
    s=t;
    q2=h_5_ceshiair5(d2,s,n);%调用积分函数五
    q3=q2-q1;
    q4=q4+q3;    
end 
if d3==1 && t0>=(period1-inter1)
    s=(d2+1/2)*T-d;
    q1=h_3_ceshiair3(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_3_ceshiair3(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T;
    q1=h_5_ceshiair5(d2,s,n);%调用积分函数五
    s=t;
    q2=h_5_ceshiair5(d2,s,n);%调用积分函数五
    q3=q2-q1;
    q4=q4+q3;    
    s=(d2+1)*T-d;
    q1=h_4_ceshiair4(d2,s,n);%调用积分函数四
    s=t;
    q2=h_4_ceshiair4(d2,s,n);%调用积分函数四
    q3=q2-q1;
    q4=q4+q3;
end 

fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q5=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/D_d1/(n*pi/L)^2;
q6=eta1*2*(-1)^n/D_d1/(n*pi/L)^3*t;
q7=eta1*2*(-1)^n/D_d1^2/(n*pi/L)^5*(exp(-D_d1*(n*pi/L)^2*t)-1);
q=q4+q5+q6+q7;
end
%第三部分
function q = IL4_3(x,t)
global A_0
global L
global B_2  %IL4的分泌速率 
global period1  %period1为一个周期
global inter1
global m
d_0=m;%初始右边界药物的条件
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q0=0;
q4=0;
for k=0:1:d2-1
    s=k*T;
    q1=h_7_ceshiair7(s);%调用积分函数七
    s=(k+1/2)*T-d;
    q2=h_7_ceshiair7(s);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;
    s=(k+1/2)*T-d;
    q1=h_8_ceshiair8(k,s);%调用积分函数八
    s=(k+1/2)*T;
    q2=h_8_ceshiair8(k,s);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3;
    s=(k+1)*T-d;
    q1=h_9_ceshiair9(k,s);%调用积分函数九
    s=(k+1)*T;
    q2=h_9_ceshiair9(k,s);%调用积分函数九
    q3=q2-q1;
    q4=q4+q3;    
end
%判断t最后走到第几个周期积分对应着对应d(t)的最后一部分的积分
if d3==0 && t0<=(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);%调用积分函数七
    s=t;
    q2=h_7_ceshiair7(s);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==0 && t0>(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);%调用积分函数七
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s);%调用积分函数八
    s=t;
    q2=h_8_ceshiair8(d2,s);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3;    
end 
if d3==1 && t0<=(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);%调用积分函数七
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s);%调用积分函数八
    s=(d2+1/2)*T;
    q2=h_8_ceshiair8(d2,s);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3; 
end 
if d3==1 && t0>(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);%调用积分函数七
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s);%调用积分函数七
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s);%调用积分函数八
    s=(d2+1/2)*T;
    q2=h_8_ceshiair8(d2,s);%调用积分函数八
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1)*T-d;
    q1=h_9_ceshiair9(d2,s);%调用积分函数九
    s=t;
    q2=h_9_ceshiair9(d2,s);%调用积分函数九
    q3=q2-q1;
    q4=q4+q3; 
end 
for n=1:1:20
    q0=q0+TnIL4_3(t,n)*sin(n*pi*x/L);
end
q=A_0+B_2/x*q0+B_2*d_0*t+B_2*q4;
end
