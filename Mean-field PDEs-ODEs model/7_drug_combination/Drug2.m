function aaa = Drug2(x,t,select,dmax,dmax_1)
% global eta_d    %药物的降解率
% global D_d    %药物的扩散率
% global delt_t  %模拟时候的时间间隔
% global L      %模拟时候横轴x的大小
% global d_max         %药物边界值，等于3
% global period        %这个是用于边界用药时候有一个周期，例如7天一个周期
% global inter         %阶梯函数的换化间隔，例如用一天把药物量为0渐变成药物量为6，一天就是这个inter
% global period1       %真实的的时间周期，由period×delt_t得到
% global inter1        %真实的的时间间隔，由inter×delt_t得到
% global TimeLength    %总的时间
% global m             %边界条件三设计的药物的量
% Parameters_nondimensional_4_gai_jin_fangcheng_4
% delt_t=1/91.5926/24;      %模拟时候的时间间隔
% TimeLength=200;
% m=3;%这个是第三种情况下，阶梯函数的最高值，代表喂药一段时间的量为m，一段时间的量为0
% period=20;           %小周期，例如20天用药量6，20天用药量0
% period1=period*delt_t*24;  %真正周期
% inter=1;             %间隔
% inter1=inter*delt_t*24;   %真正间隔
% L=100;         %模拟时候坐标轴的总长度
% eta_d=eta_d;
% D_d=D_d;%药物的扩散率
% d_max=3;%药物边界的量
global d_max m0 b1 b2

d_max=dmax;
m0=dmax;

if x==0
    x=0.0000001;
end
switch select
    case 0
        aaa=0;
    case 1
        
        aaa=Drug_1(x,t);
    case 2
        aaa=Drug_2(x,t);
    case 3
        aaa=Drug_3(x,t);
    case 4
        aaa=Drug_4(x,t);
    case 5
        aaa=Drug_5(x,t) ;
end
end
function q = Drug_1(x,t)%常数给药
global d_max
global eta_d
global D_d
global L   %x的最大值
C_d=d_max;
drugt1=0;
for n=1:1:100
    fi_n=(2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n);
    Tn1=(fi_n*exp(-(D_d*(n*pi/L)^2+eta_d)*t)+2*eta_d*C_d*L/(n*pi)*(-1)^n/(D_d*(n*pi/L)^2+eta_d)*(1-exp(-(D_d*(n*pi/L)^2+eta_d)*t)));
    Tn11=Tn1*sin(n*pi*x/L);
    drugt1=drugt1+Tn11;
end
q=1/x*drugt1+C_d; %药物的解析解
end

function q = Drug_2(x,t)
global d_max
global eta_d
global D_d
global L   %x的最大值
global  period_2
% d_initial=d_max/L;
w=2*pi/period_2;
% w=pi/(TimeLength*delt_t/(period1/delt_t));
C_d=d_max;
d22=C_d/2*(sin(w*t)+2);%d(t)
drugt2=0;
for n=1:1:100
    fi_n=(2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n);
    Tn2_1=fi_n*exp(-(D_d*(n*pi/L)^2+eta_d)*t);
    Tn2_2=2*eta_d*C_d*L/(n*pi)*(-1)^n/(D_d*(n*pi/L)^2+eta_d)*(1-exp(-(D_d*(n*pi/L)^2+eta_d)*t));
    Tn2_3=eta_d*C_d*(D_d*(n*pi/L)^2+eta_d)*L/(n*pi)*(-1)^n/((D_d*(n*pi/L)^2+eta_d)^2+w^2)*sin(w*t);
    Tn2_4=-eta_d*C_d*w*L/(n*pi)*(-1)^n/((D_d*(n*pi/L)^2+eta_d)^2+w^2)*cos(w*t);
    Tn2_5=eta_d*C_d*w*L/(n*pi)*(-1)^n/((D_d*(n*pi/L)^2+eta_d)^2+w^2)*exp(-(D_d*(n*pi/L)^2+eta_d)*t);
    Tn2_6=C_d*w*(D_d*(n*pi/L)^2+eta_d)*L/(n*pi)*(-1)^n/((D_d*(n*pi/L)^2+eta_d)^2+w^2)*cos(w*t);
    Tn2_7=-C_d*w*(D_d*(n*pi/L)^2+eta_d)*L/(n*pi)*(-1)^n/((D_d*(n*pi/L)^2+eta_d)^2+w^2)*exp(-(D_d*(n*pi/L)^2+eta_d)*t);
    Tn2_8=C_d*w^2*L/(n*pi)*(-1)^n/((D_d*(n*pi/L)^2+eta_d)^2+w^2)*sin(w*t);
    Tn2=Tn2_1+Tn2_2+Tn2_3+Tn2_4+Tn2_5+Tn2_6+Tn2_7+Tn2_8;
    drugt2=drugt2+Tn2*sin(n*pi*x/L);
end
q=1/x*drugt2+d22;
end


%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分
%积分函数一，对应h1(t,t1)  t1是to
function q=h_1_ceshiair2_1(k,t1,t,n) %g1k(t)
global D_d
global L
global inter1
global period1
global m0
global eta_d
T=2*period1;
q1=m0/(D_d*(n*pi/L)^2+eta_d)*(exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-k*T)));
q=2*eta_d*L*(-1)^n/n/pi*q1;
end
%积分函数二，对应h2(t,t1)
function q=h_2_ceshiair2_2(k,t1,t,n)%h2d2(t)和g2k(t)的统一格式
global D_d
global L
global inter1
global period1
global m0
global eta_d
T=2*period1;
d=inter1;
q11=m0/d^2*(1+2/d*(t1-(k+1/2)*T+d))*(t1-(k+1/2)*T)^2*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))-m0*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-(k+1/2)*T+d));
q1=q11/(D_d*(n*pi/L)^2+eta_d);
q21=(2*m0/d^3*(t1-(k+1/2)*T)^2+2*m0/d^2*(t1-(k+1/2)*T)*(1+2/d*(t1-(k+1/2)*T+d)))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1));
q2=-q21/(D_d*(n*pi/L)^2+eta_d)^2;
q31=(8*m0/d^3*(t1-(k+1/2)*T)+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d)))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))+6*m0/d^2*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-(k+1/2)*T+d));
q3=q31/(D_d*(n*pi/L)^2+eta_d)^3;
q41=12*m0/d^3*(exp(-(D_d*(n*pi/L)^2+eta_d)*(t-(k+1/2)*T+d))-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1)));
q4=q41/(D_d*(n*pi/L)^2+eta_d)^4;
q=(q1+q2+q3+q4)*2*eta_d*L*(-1)^n/n/pi;
end
%积分函数三，对应h2(t,t1)
function q=h_3_ceshiair2_3(k,t1,t,n)
global D_d
global L
global inter1
global period1
global m0
global eta_d
T=2*period1;
d=inter1;
q11=m0/d^2*(1-2/d*(t1-(k+1)*T))*(t1-(k+1)*T+d)^2*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1));
q1=q11/(D_d*(n*pi/L)^2+eta_d);
q21=-(-2*m0/d^3*(t1-(k+1)*T+d)^2+2*m0/d^2*(t1-(k+1)*T+d)*(1-2/d*(t1-(k+1)*T)))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1));
q2=q21/(D_d*(n*pi/L)^2+eta_d)^2;
q31=(-8*m0/d^3*(t1-(k+1)*T+d)+2*m0/d^2*(1-2/d*(t1-(k+1)*T)))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))-6*m0/d^2*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-(k+1)*T+d));
q3=q31/(D_d*(n*pi/L)^2+eta_d)^3;
q41=12*m0/d^3*(exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-(k+1)*T+d)));
q4=q41/(D_d*(n*pi/L)^2+eta_d)^4;
q=(q1+q2+q3+q4)*2*eta_d*L*(-1)^n/n/pi;
end
%积分函数四，对应h4(t,t1)
function q=h_4_ceshiair2_4(k,t1,t,n)
global D_d
global L
global inter1
global period1
global m0
global eta_d
T=2*period1;
d=inter1;
q11=(2*m0/d^3*(t1-(k+1/2)*T)^2+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d))*(t1-(k+1/2)*T))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1));
q1=q11/(D_d*(n*pi/L)^2+eta_d);
q21=-((8*m0/d^3*(t1-(k+1/2)*T)+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d)))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))+6*m0/d^2*exp(-(D_d*(n*pi/L)^2+eta_d))*(t-(k+1/2)*T+d));
q2=q21/(D_d*(n*pi/L)^2+eta_d)^2;
q31=12*m0/d^3*(exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-(k+1/2)*T+d)));
q3=q31/(D_d*(n*pi/L)^2+eta_d)^3;
q=(q1+q2+q3)*2*L*(-1)^n/n/pi;
end
%积分函数五，对应h5(t,t1)
function q=h_5_ceshiair2_5(k,t1,t,n)
global D_d
global L
global inter1
global period1
global m0
global eta_d
T=2*period1;
d=inter1;
q11=(-2*m0/d^3*(t1-(k+1)*T+d)^2+2*m0/d^2*(1-2/d*(t1-(k+1)*T))*(t1-(k+1)*T+d))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1));
q1=q11/(D_d*(n*pi/L)^2+eta_d);
q21=(8*m0/d^3*(t1-(k+1)*T+d)-2*m0/d^2*(1-2/d*(t1-(k+1)*T)))*exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))+6*m0/d^2*exp(-(D_d*(n*pi/L)^2+eta_d))*(t-(k+1)*T+d);
q2=q21/(D_d*(n*pi/L)^2+eta_d)^2;
q31=-12*m0/d^3*(exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t1))-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-(k+1)*T+d)));
q3=q31/(D_d*(n*pi/L)^2+eta_d)^3;
q=(q1+q2+q3)*2*L*(-1)^n/n/pi;
end
%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q = Tn34(t,n)
global D_d;
global eta_d;
global L
global period1
global inter1
global m0
T=2*period1;
d=inter1;
q=0;
q1=0;
q2=0;
q3=0;
q4=0;
d_0=m0;%初始条件为d_initial*r，应该要除以L
d_initial=d_0/L;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
if d2~=0 %避免0算进去了
    for i=0:1:(d2-1)
        t1=i*T;
        q1=h_1_ceshiair2_1(i,t1,t,n); %0
        t1=(i+1/2)*T-d;
        q2=h_1_ceshiair2_1(i,t1,t,n);
        q3=q2-q1;
        q4=q4+q3; %g1i(t)               
        
        t1=(i+1/2)*T-d;
        q1=h_2_ceshiair2_2(i,t1,t,n);       
        t1=(i+1/2)*T;
        q2=h_2_ceshiair2_2(i,t1,t,n);
        q3=q2-q1;
        q4=q4+q3;
        
        t1=(i+1)*T-d;
        q1=h_3_ceshiair2_3(i,t1,t,n);       
        t1=(i+1)*T;
        q2=h_3_ceshiair2_3(i,t1,t,n);
        q3=q2-q1;
        q4=q4+q3;
        
        t1=(i+1/2)*T-d;
        q1=h_4_ceshiair2_4(i,t1,t,n);       
        t1=(i+1/2)*T;
        q2=h_4_ceshiair2_4(i,t1,t,n);
        q3=q2-q1;
        q4=q4+q3;
        
        t1=(i+1)*T-d;
        q1=h_5_ceshiair2_5(i,t1,t,n);   
        t1=(i+1)*T;
        q2=h_5_ceshiair2_5(i,t1,t,n);
        q3=q2-q1;
        q4=q4+q3;
    end
end
if d3==0 && t0<(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=t;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end
if d3==0 && t0>=(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=(d2+1/2)*T-d;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1/2)*T-d;
    q1=h_2_ceshiair2_2(d2,t1,t,n);
    t1=t;
    q2=h_2_ceshiair2_2(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1/2)*T-d;
    q1=h_4_ceshiair2_4(d2,t1,t,n);
    t1=t;
    q2=h_4_ceshiair2_4(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end
if d3==1 && t0<(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=(d2+1/2)*T-d;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1/2)*T-d;
    q1=h_2_ceshiair2_2(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_2_ceshiair2_2(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1/2)*T-d;
    q1=h_4_ceshiair2_4(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_4_ceshiair2_4(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end
if d3==1 && t0>=(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=(d2+1/2)*T-d;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1/2)*T-d;
    q1=h_2_ceshiair2_2(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_2_ceshiair2_2(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1)*T-d;
    q1=h_3_ceshiair2_3(d2,t1,t,n);
    t1=t;
    q2=h_3_ceshiair2_3(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1/2)*T-d;
    q1=h_4_ceshiair2_4(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_4_ceshiair2_4(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1)*T-d;
    q1=h_5_ceshiair2_5(d2,t1,t,n);
    t1=t;
    q2=h_5_ceshiair2_5(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end
fi_n=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n);
q=fi_n*exp(-(D_d*(n*pi/L)^2+eta_d)*t)+q4;
end


%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q = Drug_3(x,t)
global L
global period1
global inter1
global m0
q0=0;
d1=floor(t/period1);%floor向负无穷取整  T=2*period1
d2=floor(d1/2);
d3=rem(d1,2);%rem取余数
t0=t-d1*period1;
xk=(d1+1)*period1-inter1;%inter1就是d
xk1=(d1+1)*period1;
if d3==0 && t0<=(period1-inter1)
    dnt=m0;
end
if d3==0 && t0>(period1-inter1)
    yk=m0;
    %yk1=0;
    dnt=(1+2*(t-xk)/(xk1-xk))*((t-xk1)/(xk-xk1))^2*yk+0+0+0;
end
if d3==1 && t0<=(period1-inter1)
    dnt=0;
end
if d3==1 && t0>(period1-inter1)
    %yk=0;
    yk1=m0;
    dnt=0+(1+2*(t-xk1)/(xk-xk1))*((t-xk)/(xk1-xk))^2*yk1+0+0;
end
for n=1:1:20
    q0=q0+Tn34(t,n)*sin(n*pi*x/L);
end
q=q0/x+dnt;
end


%下面是用药停药函数的函数，药物情况4
function q=h_0_air(t,n)
global D_d
global L
global inter1
global eta_d
global b1
d=inter1;
q11=b1/(D_d*(n*pi/L)^2+eta_d);
q12=1-exp(-(D_d*(n*pi/L)^2+eta_d)*t);
q=q11*q12;
end

function q=h_1_air(t52,t53,t,n)
global D_d
global L
global inter1
global eta_d
global b1 b2
d=inter1;
q01=b1/(D_d*(n*pi/L)^2+eta_d);
q02=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))-exp(-(D_d*(n*pi/L)^2+eta_d)*t);
q0=q01*q02;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*b1/d^2*(t52-t53)^2;
q13=(b1/d^2*(1+2/d*(t-t52))*(t-t53)^2+b2/d^2*(1-2/d*(t-t53))*(t-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q23=(2*b1/d^3*(t-t53)^2+2*b1/d^2*(1+2/d*(t-t52))*(t-t53)-2*b2/d^3*(t-t52)^2+2*b2/d^2*(1-2/d*(t-t53))*(t-t52));
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=(8*b1/d^3*(t-t53)+2*b1/d^2*(1+2/d*(t-t52))-8*b2/d^3*(t-t52)+2*b2/d^2*(1-2/d*(t-t53)));
q3=q31*(q33-q32);
q41=-1/(D_d*(n*pi/L)^2+eta_d)^4;
q42=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q43=(12*b1/d^3-12*b2/d^3);
q4=q41*(q43-q42);
q=q0+q1+q2+q3+q4;
end

function q=h_2_air(t52,t53,t,n)
global D_d
global L
global inter1
global eta_d
global b1 b2
d=inter1;
q01=b1/(D_d*(n*pi/L)^2+eta_d);
q02=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))-exp(-(D_d*(n*pi/L)^2+eta_d)*t);
q0=q01*q02;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*b1/d^2*(t52-t53)^2;
q13=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(b2/d^2*(t53-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q23=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q3=q31*(q33-q32);
q41=-1/(D_d*(n*pi/L)^2+eta_d)^4;
q42=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q43=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(12*b1/d^3-12*b2/d^3);
q4=q41*(q43-q42);
q51=b2/(D_d*(n*pi/L)^2+eta_d);
q52=1-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53));
q5=q51*q52;
q=q0+q1+q2+q3+q4+q5;
end


function q=h_3_air(t52,t53,t,n)
global D_d
global L
global inter1
global period1
global eta_d
global b1 b2
d=inter1;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q13=(2*b1/d^3*(t-t53)^2+2*b1/d^2*(1+2/d*(t-t52))*(t-t53)-2*b2/d^3*(t-t52)^2+2*b2/d^2*(1-2/d*(t-t53))*(t-t52));
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=(8*b1/d^3*(t-t53)+2*b1/d^2*(1+2/d*(t-t52))-8*b2/d^3*(t-t52)+2*b2/d^2*(1-2/d*(t-t53)));
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q33=(12*b1/d^3-12*b2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end
%积分函数二，对应h4
function q=h_4_air(t52,t53,t,n)
global D_d
global L
global inter1
global period1
global eta_d
global b1 b2
d=inter1;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q13=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q33=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(12*b1/d^3-12*b2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end
%%这个是用药为阶梯函数，一段时间用药b1，一段时间用药b2，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q = Tn4(t,n)
global D_d;
global eta_d;
global L
global period1
global inter1
global b1 b2
global t51 t52 t53
d_0=b1;%初始条件为d_initial*r，应该要除以L
d_initial=d_0/L;
if t>=t51 && t<t52
    q0=h_0_air(t,n)*eta_d*2*L/(n*pi)*(-1)^n;
    q1=0;
elseif t>=t52 && t<t53
    q0=h_1_air(t52,t53,t,n)*eta_d*2*L/(n*pi)*(-1)^n;
    q1=h_3_air(t52,t53,t,n)*2*L/(n*pi)*(-1)^n;
elseif t>=t53
    q0=h_2_air(t52,t53,t,n)*eta_d*2*L/(n*pi)*(-1)^n;
    q1=h_4_air(t52,t53,t,n)*2*L/(n*pi)*(-1)^n;
else
    msg='Error occured';
    error(msg)
end
q2=q0+q1;
q=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n)*exp(-(D_d*(n*pi/L)^2+eta_d)*t)+q2;
end


%%这个是用药为阶梯l两段函数，一段时间用药b1，一段时间用药b2，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q = Drug_4(x,t)
global L
global t51 t52 t53
global inter1
global b1 b2
d=inter1;
if t>=t51 && t<t52
    dnt=b1;
elseif t>=t52 && t<t53
    dnt=b1/d^2*(1+2/d*(t-t52))*(t-t53)^2+b2/d^2*(1-2/d*(t-t53))*(t-t52)^2;
elseif t>=t53
    dnt=b2;
else
    msg='Error occured';
    error(msg)
end
q0=0;
for n=1:1:20
    q0=q0+Tn4(t,n)*sin(n*pi*x/L);
end
q=q0/x+dnt;
end



%下面是用药停药函数的函数，药物情况5
function q=p1i_air(t51,t,n,bi)
global D_d
global L
global inter1
global eta_d
d=inter1;
q11=bi/(D_d*(n*pi/L)^2+eta_d);
q12=1-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t51));
q=q11*q12;
end
%积分函数二，对应h4
function q=p2i_air(t52,t53,t,n,bi1,bi2)
global D_d
global L
global inter1
global period1
global eta_d
d=inter1;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*bi1/d^2*(t52-t53)^2;
q13=(bi1/d^2*(1+2/d*(t-t52))*(t-t53)^2+bi2/d^2*(1-2/d*(t-t53))*(t-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q23=(2*bi1/d^3*(t-t53)^2+2*bi1/d^2*(1+2/d*(t-t52))*(t-t53)-2*bi2/d^3*(t-t52)^2+2*bi2/d^2*(1-2/d*(t-t53))*(t-t52));
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=(8*bi1/d^3*(t-t53)+2*bi1/d^2*(1+2/d*(t-t52))-8*bi2/d^3*(t-t52)+2*bi2/d^2*(1-2/d*(t-t53)));
q3=q31*(q33-q32);
q41=-1/(D_d*(n*pi/L)^2+eta_d)^4;
q42=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=(12*bi1/d^3-12*bi2/d^3);
q4=q41*(q43-q42);
q=q1+q2+q3+q4;
end
function q=p3i_air(t52,t53,t,n,bi1,bi2)
global D_d
global L
global inter1
global period1
global eta_d
d=inter1;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q13=(2*bi1/d^3*(t-t53)^2+2*bi1/d^2*(1+2/d*(t-t52))*(t-t53)-2*bi2/d^3*(t-t52)^2+2*bi2/d^2*(1-2/d*(t-t53))*(t-t52));
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=(8*bi1/d^3*(t-t53)+2*bi1/d^2*(1+2/d*(t-t52))-8*bi2/d^3*(t-t52)+2*bi2/d^2*(1-2/d*(t-t53)));
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=(12*bi1/d^3-12*bi2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end
function q=q1i_air(t51,t52,t,n,bi)
global D_d
global L
global inter1
global eta_d
d=inter1;
q11=bi/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))-exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t51));
q=q11*q12;
end
%积分函数二，对应h4
function q=q2i_air(t52,t53,t,n,bi1,bi2)
global D_d
global L
global inter1
global period1
global eta_d
d=inter1;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*bi1/d^2*(t52-t53)^2;
q13=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(bi2/d^2*(t53-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q23=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q3=q31*(q33-q32);
q41=-1/(D_d*(n*pi/L)^2+eta_d)^4;
q42=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(12*bi1/d^3-12*bi2/d^3);
q4=q41*(q43-q42);
q=q1+q2+q3+q4;
end
function q=q3i_air(t52,t53,t,n,bi1,bi2)
global D_d
global L
global inter1
global period1
global eta_d
d=inter1;
q11=1/(D_d*(n*pi/L)^2+eta_d);
q12=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q13=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q1=q11*(q13-q12);
q21=-1/(D_d*(n*pi/L)^2+eta_d)^2;
q22=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q2=q21*(q23-q22);
q31=1/(D_d*(n*pi/L)^2+eta_d)^3;
q32=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=exp(-(D_d*(n*pi/L)^2+eta_d)*(t-t53))*(12*bi1/d^3-12*bi2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end
%%这个是用药为阶梯函数，一段时间用药b1，一段时间用药b2，一段时间用药b3,一段时间用药b4,连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q = Tn5(t,n)
global D_d;
global eta_d;
global L
global period1
global inter1
global b1 b2 b3 b4
global t51 t52 t53 t54 t55 t56 t57 t58
global b t5
d_0=b(1);%初始条件为d_initial*r，应该要除以L
d_initial=d_0/L;


aa1=find(t>=t5);
aa2=aa1(end);
aa3=floor(aa2/2);
aa4=aa2-aa3*2;
q11=0;q12=0;q13=0;
if aa3==0
    q1=eta_d*p1i_air(t5(1),t,n,b(1));
elseif aa3==1
    if aa4==0
        q11=q11+eta_d*(q1i_air(t5(1),t5(2),t,n,b(1)));
        q1=q11+eta_d*p2i_air(t5(2),t5(3),t,n,b(1),b(2))+p3i_air(t5(2),t5(3),t,n,b(1),b(2));
    else
        q11=eta_d*(q1i_air(t5(1),t5(2),t,n,b(1)));
        q12=eta_d*q2i_air(t5(2),t5(3),t,n,b(1),b(2));
        q13=q3i_air(t5(2),t5(3),t,n,b(1),b(2));
        q1=q11+q12+q13+eta_d*p1i_air(t5(3),t,n,b(2));
    end
else
    if aa4==0
        q11=q11+eta_d*(q1i_air(t5(1),t5(2),t,n,b(1)));
        for i=1:1:aa3-1
            q12=q12+eta_d*q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
            q13=q13+q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
            q11=q11+eta_d*(q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1)));
        end
        q1=q11+q12+q13+eta_d*p2i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1))+p3i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1));
    else
        q11=q11+eta_d*(q1i_air(t5(1),t5(2),t,n,b(1)));
        for i=1:1:aa3-1
            q12=q12+eta_d*q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
            q13=q13+q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
            q11=q11+eta_d*(q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1)));
        end
        q12=q12+eta_d*q2i_air(t5(aa2-1),t5(aa2),t,n,b(aa3),b(aa3+1));
        q13=q13+q3i_air(t5(aa2-1),t5(aa2),t,n,b(aa3),b(aa3+1));
        q1=q11+q12+q13+eta_d*p1i_air(t5(aa2),t,n,b(aa3+1));
    end
end








% if t>=t51 && t<t52
%     q1=eta_d*p1i_air(t51,t,n,b1);
% elseif t>=t52 && t<t53
%     q11=eta_d*q1i_air(t51,t52,t,n,b1);
%     q1=q11+eta_d*p2i_air(t52,t53,t,n,b1,b2)+p3i_air(t52,t53,t,n,b1,b2);
% elseif t>=t53 && t<t54
%     q11=eta_d*q1i_air(t51,t52,t,n,b1);
%     q12=eta_d*q2i_air(t52,t53,t,n,b1,b2);
%     q13=q3i_air(t52,t53,t,n,b1,b2);
%     q1=q11+q12+q13+eta_d*p1i_air(t53,t,n,b2);
% elseif t>=t54 && t<t55
%     q11=eta_d*(q1i_air(t51,t52,t,n,b1)+q1i_air(t53,t54,t,n,b2));
%     q12=eta_d*q2i_air(t52,t53,t,n,b1,b2);
%     q13=q3i_air(t52,t53,t,n,b1,b2);
%     q1=q11+q12+q13+eta_d*p2i_air(t54,t55,t,n,b2,b3)+p3i_air(t54,t55,t,n,b2,b3);
% elseif t>=t55 && t<t56
%     q11=eta_d*(q1i_air(t51,t52,t,n,b1)+q1i_air(t53,t54,t,n,b2));
%     q12=eta_d*(q2i_air(t52,t53,t,n,b1,b2)+q2i_air(t54,t55,t,n,b2,b3));
%     q13=(q3i_air(t52,t53,t,n,b1,b2)+q3i_air(t54,t55,t,n,b2,b3));
%     q1=q11+q12+q13+eta_d*p1i_air(t55,t,n,b3);
% elseif t>=t56 && t<t57
%     q11=eta_d*(q1i_air(t51,t52,t,n,b1)+q1i_air(t53,t54,t,n,b2)+q1i_air(t55,t56,t,n,b3));
%     q12=eta_d*(q2i_air(t52,t53,t,n,b1,b2)+q2i_air(t54,t55,t,n,b2,b3));
%     q13=(q3i_air(t52,t53,t,n,b1,b2)+q3i_air(t54,t55,t,n,b2,b3));
%     q1=q11+q12+q13+eta_d*p2i_air(t56,t57,t,n,b3,b4)+p3i_air(t56,t57,t,n,b3,b4);
% elseif t>=t57 && t<=t58
%     q11=eta_d*(q1i_air(t51,t52,t,n,b1)+q1i_air(t53,t54,t,n,b2)+q1i_air(t55,t56,t,n,b3));
%     q12=eta_d*(q2i_air(t52,t53,t,n,b1,b2)+q2i_air(t54,t55,t,n,b2,b3)+q2i_air(t56,t57,t,n,b3,b4));
%     q13=(q3i_air(t52,t53,t,n,b1,b2)+q3i_air(t54,t55,t,n,b2,b3)+q3i_air(t56,t57,t,n,b3,b4));
%     q1=q11+q12+q13+eta_d*p1i_air(t57,t,n,b4);
% else
%     msg='Error occured';
%     error(msg)
% end

fi_n=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n);
q4=2*L/(n*pi)*(-1)^n*q1;
q=fi_n*exp(-(D_d*(n*pi/L)^2+eta_d)*t)+q4;
end


%%这个是用药为阶梯l两段函数，一段时间用药b1，一段时间用药b2，一段时间用药b3,一段时间用药b4，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q = Drug_5(x,t)
global L
global t51 t52 t53 t54 t55 t56 t57 t58
global inter1
global b1 b2 b3 b4
global b t5
d=inter1;
aa1=find(t>=t5);
aa2=aa1(end);
aa3=floor(aa2/2);%朝负无穷大方向取整
aa4=aa2-aa3*2;

if aa3==0
    dnt=b(1);
else
    if aa4==0
        dnt=b(aa3)/d^2*(1+2/d*(t-t5(aa2)))*(t-t5(aa2+1))^2+b(aa3+1)/d^2*(1-2/d*(t-t5(aa2+1)))*(t-t5(aa2))^2;
    else
        dnt=b(aa3+1);
    end
end







% if t>=t51 && t<t52
%     dnt=b1;
% elseif t>=t52 && t<t53
%     dnt=b1/d^2*(1+2/d*(t-t52))*(t-t53)^2+b2/d^2*(1-2/d*(t-t53))*(t-t52)^2;
% elseif t>=t53 && t<t54
%     dnt=b2;
% elseif t>=t54 && t<t55
%     dnt=b2/d^2*(1+2/d*(t-t54))*(t-t55)^2+b3/d^2*(1-2/d*(t-t55))*(t-t54)^2;
% elseif t>=t55 && t<t56
%     dnt=b3;
% elseif t>=t56 && t<t57
%     dnt=b3/d^2*(1+2/d*(t-t56))*(t-t57)^2+b4/d^2*(1-2/d*(t-t57))*(t-t56)^2;
% elseif t>=t57 && t<=t58
%     dnt=b4;
% else
%     msg='Error occured';
%     error(msg)
% end

q0=0;
for n=1:1:20
    q0=q0+Tn5(t,n)*sin(n*pi*x/L);
end
q=q0/x+dnt;
end