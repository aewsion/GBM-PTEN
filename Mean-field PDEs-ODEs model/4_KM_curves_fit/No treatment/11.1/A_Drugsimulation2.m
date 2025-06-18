function aaa=A_Drugsimulation2(x,t,select,dmax,dmax_1)
% global eta1    %药物的降解率
% global D_d1    %药物的扩散率
% global delt_t  %模拟时候的时间间隔
% global L      %模拟时候横轴x的大小
% global B_3  %这个是需要估计的值,是IL4分泌速率相关的
% global A_0   %这个是需要估计的值,是耐药细胞因子初始值
% global d_max         %药物边界值，等于3
% global period        %这个是用于边界用药时候有一个周期，例如7天一个周期
% global inter         %阶梯函数的换化间隔，例如用一天把药物量为0渐变成药物量为6，一天就是这个inter
% global period1       %真实的的时间周期，由period×delt_t得到
% global inter1        %真实的的时间间隔，由inter×delt_t得到
% global TimeLength    %总的时间
% global m             %边界条件三设计的药物的量
% Parameters_nondimensional_4_gai_jin_fangcheng_4
% TimeLength=200;
% m=3;%这个是第三种情况下，阶梯函数的最高值，代表喂药一段时间的量为m，一段时间的量为0
% delt_t=1/91.5926/24;      %模拟时候的时间间隔
% period=20;           %小周期，例如20天用药量6，20天用药量0
% period1=period*delt_t*24;  %真正周期
% inter=1;             %间隔
% inter1=inter*delt_t*24;   %真正间隔
% L=100;         %模拟时候坐标轴的总长度
% eta1=eta_d;
% D_d1=D_d;%药物的扩散率
% d_max=3;%药物边界的量
global d_max m0 b1 b2
% if isempty(dmax)
% else
    d_max=dmax;
    m0=dmax;
%     b1=dmax;
% end
% if nargin >=5
%     b2=dmax_1;
% end
if x==0
    x=0.0000001;
end
switch select
    case 0
        aaa=0;
    case 1
        aaa=A_1(x,t);
    case 2
        aaa=A_2(x,t);
    case 3
        aaa=A_3(x,t);
    case 4
        aaa=A_4(x,t);
    case 5
        aaa=A_5(x,t);
end
end
%下面是计算在位置x和时间t时由药物CSF1R_I产生的药物产生的A的浓度

%%药物边界条件为常量3的时候产生的A的量
%第一部分
function q = TnA_1(t,n)
global d_max
global D_d1;
global eta1;
global L
C_d=d_max;
% d_initial=d_max/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
fi_n=2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1)+2*eta1*C_d*(L/n/pi)*(-1)^n/(D_d1*(n*pi/L)^2+eta1)*(t+(exp(-(D_d1*(n*pi/L)^2+eta1)*t)-1)/(D_d1*(n*pi/L)^2+eta1));
end
%第二部分
function q = A_1(x,t)
global d_max
global A_0
global L
global B_3  %A的分泌速率 
q0=0;
C_d=d_max;
for n=1:1:100
    q0=q0+sin(n*pi/L*x)*TnA_1(t,n);
end
q=A_0+B_3/x*q0+B_3*C_d*t;
end

%%药物边界条件为d22=d_max/2*(sin(w*t)+2)=3/2*(sin(w*t)+2)的时候产生的A的量
%第一部分
function q = TnA_2(t,n)
global d_max
global D_d1;
global eta1;
global L
global period_2
% w=pi/(TimeLength*delt_t/20);
w=2*pi/period_2;
C_d=d_max;
% d_initial=d_max/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
fi_n=(2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n);
q1=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q2=2*eta1*C_d*(L/n/pi)*(-1)^n/(D_d1*(n*pi/L)^2+eta1)*(t+(exp(-(D_d1*(n*pi/L)^2+eta1)*t)-1)/(D_d1*(n*pi/L)^2+eta1));
q3=eta1*C_d*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-cos(w*t))/w;
q4=-eta1*C_d*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*sin(w*t);
q5=eta1*C_d*w*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q6=C_d*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*sin(w*t);
q7=-C_d*w*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q8=C_d*w*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-cos(w*t));
q=q1+q2+q3+q4+q5+q6+q7+q8;
end

%第二部分
function q = A_2(x,t)
global d_max
global A_0
global L
global B_3  %A的分泌速率 
global period_2;
C_d=d_max;
q0=0;

w=2*pi/period_2;
for n=1:1:100
    q0=q0+sin(n*pi*x/L)*TnA_2(t,n);
end
q=A_0+B_3/x*q0+B_3*C_d/2/w*(1-cos(w*t))+B_3*C_d*t;
end

%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分


%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数三
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%第一部分
function q=h_11_ceshiair2_11(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q1=2*eta1*L*(-1)^n/(n*pi)*m0/(D_d1*(n*pi/L)^2+eta1);
q2=s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-k*T))/(D_d1*(n*pi/L)^2+eta1);
q=q1*q2;
end
function q=h_12_ceshiair2_12(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0/d^2*((s-(k+1/2)*T)^3+1/(2*d)*(s-(k+1/2)*T)^4);
q12=m0*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q1=(q11+q12)/(D_d1*(n*pi/L)^2+eta1);
q21=2*m0/(3*d^3)*(s-(k+1/2)*T)^3+2*m0/d^2*(3/2*(s-(k+1/2)*T)^2+2/(3*d)*(s-(k+1/2)*T)^3);
q2=-q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=4*m0/d^3*(s-(k+1/2)*T)^2+2*m0/d^2*(s+1/d*(s-(k+1/2)*T+d)^2);
q32=-6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q3=(q31+q32)/(D_d1*(n*pi/L)^2+eta1)^3;
q41=-12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1)+s);
q4=q41/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2+q3+q4);
end
function q=h_13_ceshiair2_13(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0/d^2*((s-(k+1)*T+d)^3-1/(2*d)*(s-(k+1)*T+d)^4);
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=-2*m0/(3*d^3)*(s-(k+1)*T+d)^3+2*m0/d^2*(3/2*(s-(k+1)*T+d)^2-2/(3*d)*(s-(k+1)*T+d)^3);
q2=-q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=-4*m0/d^3*(s-(k+1)*T+d)^2+2*m0/d^2*(s-1/d*(s-(k+1)*T)^2);
q32=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q3=(q31+q32)/(D_d1*(n*pi/L)^2+eta1)^3;
q41=12*m0/d^3*(s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1));
q4=q41/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2+q3+q4);
end
function q=h_14_ceshiair2_14(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=2*m0/(3*d^3)*(s-(k+1/2)*T)^3+2*m0/d^2*(3/2*(s-(k+1/2)*T)^2+2/(3*d)*(s-(k+1/2)*T)^3);
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=4*m0/d^3*(s-(k+1/2)*T)^2+2*m0/d^2*(s+1/d*(s-(k+1/2)*T+d)^2);
q22=-6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q2=-(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^2;
q31=12*m0/d^3*(s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1));
q3=q31/(D_d1*(n*pi/L)^2+eta1)^3;
q=q0*(q1+q2+q3);
end
function q=h_15_ceshiair2_15(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=-2*m0/(3*d^3)*(s-(k+1)*T+d)^3+2*m0/d^2*(3/2*(s-(k+1)*T+d)^2-2/(3*d)*(s-(k+1)*T+d)^3);
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=4*m0/d^3*(s-(k+1)*T+d)^2-2*m0/d^2*(s-1/d*(s-(k+1)*T)^2);
q22=-6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q2=(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^2;
q31=12*m0/d^3*(s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1));
q3=-q31/(D_d1*(n*pi/L)^2+eta1)^3;
q=q0*(q1+q2+q3);
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为对g1函数进行积分
function q=g_11_ceshiair2_11(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q1=2*eta1*L*(-1)^n/(n*pi)*m0/(D_d1*(n*pi/L)^2+eta1)^2;
q2=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-k*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q=q1*q2;
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为对g2函数进行积分
function q=g_12_ceshiair2_12(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q1=q11/(D_d1*(n*pi/L)^2+eta1)^2;
q21=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T));
q22=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q2=-(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^4;
q31=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T)));
q3=-q31/(D_d1*(n*pi/L)^2+eta1)^5;
q=q0*(q1+q2+q3);
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为对g3函数进行积分
function q=g_13_ceshiair2_13(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T));
q1=-q11/(D_d1*(n*pi/L)^2+eta1)^2;
q21=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d));
q22=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T));
q2=(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^4;
q31=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d)));
q3=-q31/(D_d1*(n*pi/L)^2+eta1)^5;
q=q0*(q1+q2+q3);
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为对g4函数进行积分
function q=g_14_ceshiair2_14(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T));
q12=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q1=(q11+q12)/(D_d1*(n*pi/L)^2+eta1)^3;
q21=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d)));
q2=-q21/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2);
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为对g5函数进行积分
function q=g_15_ceshiair2_15(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T));
q12=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d));
q1=-(q11+q12)/(D_d1*(n*pi/L)^2+eta1)^3;
q21=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d)));
q2=q21/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2);
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数七
%对应的是对d(t)阶梯函数第一部分的积分
function q1=h_7_ceshiair7(s)
global m0
q1=m0*s;
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数八
%对应的是对d(t)阶梯函数第二部分的积分
function q5=h_8_ceshiair8(k,s)
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=m0/d^2/2/d*s^4;
q2=m0/d^2/3*(3-2/d*(k+1/2)*T-4/d*(k+1/2)*T)*s^3;
q3=m0/d^2/2*(-6+6/d*(k+1/2)*T)*(k+1/2)*T*s^2;
q4=m0/d^2*(3-2/d*(k+1/2)*T)*(k+1/2)^2*T^2*s;
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数八
%对应的是对d(t)阶梯函数第二部分的积分
function q5=h_9_ceshiair9(k,s)
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=-m0/d^2/2/d*s^4;
q2=m0/d^2/3*(1+2/d*(k+1)*T-4/d*(d-(k+1)*T))*s^3;
q3=m0/d^2/2*(-2/d*(d-(k+1)*T)^2+2*(1+2/d*(k+1)*T)*(d-(k+1)*T))*s^2;
q4=m0/d^2*(1+2/d*(k+1)*T)*(d-(k+1)*T)^2*s;
q5=q1+q2+q3+q4;
end
%%药物边界条件为阶梯函数6和0的时候产生的A的量，中间用埃尔米插值连接
%第二部分
function q = TnA_3(t,n)
global D_d1;
global eta1;
global L
global period1  %period1为一个周期
global inter1
global m0
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q=0;
q4=0;
d_0=m0;%药物右边界条件，持续输入
d_initial=d_0/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%积分函数五是因为fn包含边界函数d'(t)是阶梯函数，分段积分完成了后是一个定的关于t的积分，再对t积分后和积分函数三不一样，然后才是IL的值的一部分
if d2~=0
for i=0:1:(d2-1)
    s=i*T;
    q1=h_11_ceshiair2_11(i,s,n);%调用积分函数三
    s=(i+1/2)*T-d;
    q2=h_11_ceshiair2_11(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T-d;
    q1=h_12_ceshiair2_12(i,s,n);%调用积分函数三
    s=(i+1/2)*T;
    q2=h_12_ceshiair2_12(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T-d;
    q1=h_13_ceshiair2_13(i,s,n);%调用积分函数三
    s=(i+1)*T;
    q2=h_13_ceshiair2_13(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T-d;
    q1=h_14_ceshiair2_14(i,s,n);%调用积分函数三
    s=(i+1/2)*T;
    q2=h_14_ceshiair2_14(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T-d;
    q1=h_15_ceshiair2_15(i,s,n);%调用积分函数三
    s=(i+1)*T;
    q2=h_15_ceshiair2_15(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T-d;
    q1=g_11_ceshiair2_11(i,s,n);%调用积分函数三
    s=t;
    q2=g_11_ceshiair2_11(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T;
    q1=g_12_ceshiair2_12(i,s,n);%调用积分函数三
    s=t;
    q2=g_12_ceshiair2_12(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T;
    q1=g_13_ceshiair2_13(i,s,n);%调用积分函数三
    s=t;
    q2=g_13_ceshiair2_13(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T;
    q1=g_14_ceshiair2_14(i,s,n);%调用积分函数三
    s=t;
    q2=g_14_ceshiair2_14(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T;
    q1=g_15_ceshiair2_15(i,s,n);%调用积分函数三
    s=t;
    q2=g_15_ceshiair2_15(i,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
end
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
if d3==0 && t0<(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    s=t;
    q2=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==0 && t0>=(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T-d;
    q2=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=g_11_ceshiair2_11(d2,s,n);%调用积分函数三
    s=t;
    q2=g_11_ceshiair2_11(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T-d;
    q1=h_12_ceshiair2_12(d2,s,n);%调用积分函数三
    s=t;
    q2=h_12_ceshiair2_12(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=h_14_ceshiair2_14(d2,s,n);%调用积分函数三
    s=t;
    q2=h_14_ceshiair2_14(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
end
if d3==1 && t0<(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T-d;
    q2=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T-d;
    q1=h_12_ceshiair2_12(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_12_ceshiair2_12(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=h_14_ceshiair2_14(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_14_ceshiair2_14(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=g_11_ceshiair2_11(d2,s,n);%调用积分函数三
    s=t;
    q2=g_11_ceshiair2_11(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T;
    q1=g_12_ceshiair2_12(d2,s,n);%调用积分函数三
    s=t;
    q2=g_12_ceshiair2_12(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T;
    q1=g_14_ceshiair2_14(d2,s,n);%调用积分函数三
    s=t;
    q2=g_14_ceshiair2_14(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
end 
if d3==1 && t0>=(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T-d;
    q2=h_11_ceshiair2_11(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T-d;
    q1=h_12_ceshiair2_12(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_12_ceshiair2_12(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=h_14_ceshiair2_14(d2,s,n);%调用积分函数三
    s=(d2+1/2)*T;
    q2=h_14_ceshiair2_14(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=g_11_ceshiair2_11(d2,s,n);%调用积分函数三
    s=t;
    q2=g_11_ceshiair2_11(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T;
    q1=g_12_ceshiair2_12(d2,s,n);%调用积分函数三
    s=t;
    q2=g_12_ceshiair2_12(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T;
    q1=g_14_ceshiair2_14(d2,s,n);%调用积分函数三
    s=t;
    q2=g_14_ceshiair2_14(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1)*T-d;
    q1=h_13_ceshiair2_13(d2,s,n);%调用积分函数三
    s=t;
    q2=h_13_ceshiair2_13(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1)*T-d;
    q1=h_15_ceshiair2_15(d2,s,n);%调用积分函数三
    s=t;
    q2=h_15_ceshiair2_15(d2,s,n);%调用积分函数三
    q3=q2-q1;
    q4=q4+q3;
end 

fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q5=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q=q4+q5;
end
%第三部分
function q = A_3(x,t)
global A_0
global L
global B_3  %A的分泌速率 
global period1  %period1为一个周期
global inter1
global m0
d_0=m0;%初始右边界药物的条件
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q0=0;
q4=0;
% if d2~=0
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
% end
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
    q0=q0+TnA_3(t,n)*sin(n*pi*x/L);
end
% q=A_0+B_3/x*q0+B_3*d_0*t+B_3*q4;
q=A_0+B_3/x*q0+B_3*q4;
end


%%用药情况4：这个是用药为阶梯分两段函数，一段时间用药b1，一段时间用药b2，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分

function q=h_0_airjifen(s,n)
global D_d1
global L
global inter1
global eta1
global b1
d=inter1;
q11=b1/(D_d1*(n*pi/L)^2+eta1);
q12=s+exp(-(D_d1*(n*pi/L)^2+eta1)*s)/(D_d1*(n*pi/L)^2+eta1);
q=q11*q12;
end

function q=h_1_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q01=b1/(D_d1*(n*pi/L)^2+eta1);
q02=exp(-(D_d1*(n*pi/L)^2+eta1)*s)-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52));
q03=1/(D_d1*(n*pi/L)^2+eta1);
q0=q01*q02*q03;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*b1/d^2*(t52-t53)^2;
q131=b1/d^2*(s-t53)^3/3;
q132=2*b1/d^3*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s);
q133=b2/d^2*(s-t52)^3/3;
q134=-2*b2/d^3*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s);
q1=q11*(q131+q132+q133+q134+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q231=2*b1/(3*d^3)*(s-t53)^3+b1/d^2*(s-t53)^2+4*b1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q232=-2*b2/(3*d^3)*(s-t52)^3+b2/d^2*(s-t52)^2-4*b2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q2=q21*(q231+q232+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=4*b1/d^3*(s-t53)^2+2*b1/d^2*(s+(s-t52)^2/d)-4*b2/d^3*(s-t52)^2+2*b2/d^2*(s-(s-t53)^2/d);
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q43=(12*b1/d^3-12*b2/d^3)*s;
q4=q41*(q43+q42/(D_d1*(n*pi/L)^2+eta1));
q=q0+q1+q2+q3+q4;
end

function q=h_2_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q01=b1/(D_d1*(n*pi/L)^2+eta1);
q02=exp(-(D_d1*(n*pi/L)^2+eta1)*s)-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52));
q03=1/(D_d1*(n*pi/L)^2+eta1);
q0=q01*q02*q03;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*b1/d^2*(t52-t53)^2;
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(b2/d^2*(t53-t52)^2);
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q43=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*b1/d^3-12*b2/d^3);
q44=1/(D_d1*(n*pi/L)^2+eta1);
q4=q41*(q42-q43)*q44;
q51=b2/(D_d1*(n*pi/L)^2+eta1);
q52=s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))/(D_d1*(n*pi/L)^2+eta1);
q5=q51*q52;
q=q0+q1+q2+q3+q4+q5;
end

%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为函数三的积分函数
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
function q=h_3_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q131=2*b1/(3*d^3)*(s-t53)^3+b1/d^2*(s-t53)^2+4*b1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q132=-b2/(3*d^3)*(s-t52)^3+b2*(s-t52)^2-4*b2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q1=q11*(q131+q132+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=4*b1/d^3*(s-t53)^2+2*b1/d^2*(s+(s-t52)^2/d)-4*b2/d^3*(s-t52)^2+2*b2/d^2*(s-(s-t53)^2/d);
q2=q21*(q23+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q33=(12*b1/d^3-12*b2/d^3)*s;
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q=q1+q2+q3;
end
function q=h_4_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*b1/d^3-12*b2/d^3);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q=q1+q2+q3;
end
function q=h_5_airjifen(s)
global b1
q=b1*s;
end
function q=h_6_airjifen(t52,t53,s)
global inter1
global b1 b2
d=inter1;
q1=b1/d^2*((s-t53)^3/3+2/d*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s));
q2=b2/d^2*((s-t52)^3/3-2/d*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s));
q=q1+q2;
end
function q=h_7_airjifen(s)
global b2
q=b2*s;
end
%第二部分
function q = TnA_4(t,n)
global D_d1;
global eta1;
global L
global b1
global t51 t52 t53
d_0=b1;%药物右边界条件，持续输入
d_initial=d_0/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%积分函数五是因为fn包含边界函数d'(t)是阶梯函数，分段积分完成了后是一个定的关于t的积分，再对t积分后和积分函数三不一样，然后才是IL的值的一部分
if t>=t51 && t<t52
    q4=eta1*(h_0_airjifen(t,n)-h_0_airjifen(0,n));
elseif t>=t52 && t<t53
    q41=eta1*(h_0_airjifen(t52,n)-h_0_airjifen(0,n));
    q42=eta1*(h_1_airjifen(t52,t53,t,n)-h_1_airjifen(t52,t53,t52,n));
    q43=h_3_airjifen(t52,t53,t,n)-h_3_airjifen(t52,t53,t52,n);
    q4=q41+q42+q43;
elseif t>=t53
    q41=eta1*(h_0_airjifen(t52,n)-h_0_airjifen(0,n));
    q42=eta1*(h_1_airjifen(t52,t53,t53,n)-h_1_airjifen(t52,t53,t52,n));
    q43=h_3_airjifen(t52,t53,t53,n)-h_3_airjifen(t52,t53,t52,n);
    q44=eta1*(h_2_airjifen(t52,t53,t,n)-h_2_airjifen(t52,t53,t53,n));
    q45=h_4_airjifen(t52,t53,t,n)-h_4_airjifen(t52,t53,t53,n);
    q4=q41+q42+q43+q44+q45;
else
    msg='Error occured';
    error(msg)    
end
q5=q4*2*L/n/pi*(-1)^n;
fi_n=2*d_0*L/n/pi*(-1)^(n+1)+4*d_0*L/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1)+q5;
end
%第三部分
function q = A_4(x,t)
global A_0
global L
global B_3  %A的分泌速率 
global period1  %period1为一个周期
global inter1
global b1
global t51 t52 t53
d_0=b1;%初始右边界药物的条件
d=inter1;
q0=0;
%判断t最后走到第几个周期积分对应着对应d(t)的最后一部分的积分
if t>=t51 && t<t52
    q41=h_5_airjifen(t51);
    q42=h_5_airjifen(t);
    q4=q42-q41;
elseif t>=t52 && t<t53
    q41=h_5_airjifen(t51);
    q42=h_5_airjifen(t52);
    q43=h_6_airjifen(t52,t53,t52);
    q44=h_6_airjifen(t52,t53,t);
    q4=q42-q41+q44-q43;
elseif t>=t53
    q41=h_5_airjifen(t51);
    q42=h_5_airjifen(t52);
    q43=h_6_airjifen(t52,t53,t52);
    q44=h_6_airjifen(t52,t53,t53);
    q45=h_7_airjifen(t53);
    q46=h_7_airjifen(t);
    q4=q42-q41+q44-q43+q46-q45;
else
    msg='Error occured';
    error(msg)    
end
for n=1:1:20
    q0=q0+TnA_4(t,n)*sin(n*pi*x/L);
end
% q=A_0+B_3/x*q0+B_3*d_0*t+B_3*q4;
q=A_0+B_3/x*q0+B_3*q4;
end




%下面是用药停药函数的函数，药物情况5
function q=jifen_p1i_air(t51,s,n,bi)
global D_d1
global L
global inter1
global eta1
d=inter1;
q11=bi/(D_d1*(n*pi/L)^2+eta1);
q12=s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t51))/(D_d1*(n*pi/L)^2+eta1);
q=q11*q12;
end
%积分函数二，对应h4
function q=jifen_p2i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*bi1/d^2*(t52-t53)^2;
q131=bi1/d^2*(s-t53)^3/3;
q132=2*bi1/d^3*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s);
q133=bi2/d^2*(s-t52)^3/3;
q134=-2*bi2/d^3*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s);
q1=q11*(q131+q132+q133+q134+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q231=2*bi1/(3*d^3)*(s-t53)^3+bi1/d^2*(s-t53)^2+4*bi1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q232=-2*bi2/(3*d^3)*(s-t52)^3+bi2/d^2*(s-t52)^2-4*bi2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q2=q21*(q231+q232+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=4*bi1/d^3*(s-t53)^2+2*bi1/d^2*(s+(s-t52)^2/d)-4*bi2/d^3*(s-t52)^2+2*bi2/d^2*(s-(s-t53)^2/d);
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=(12*bi1/d^3-12*bi2/d^3)*s;
q4=q41*(q43+q42/(D_d1*(n*pi/L)^2+eta1));
q=q1+q2+q3+q4;
end
function q=jifen_p3i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q131=2*bi1/(3*d^3)*(s-t53)^3+bi1/d^2*(s-t53)^2+4*bi1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q132=-bi2/(3*d^3)*(s-t52)^3+bi2*(s-t52)^2-4*bi2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q1=q11*(q131+q132+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=4*bi1/d^3*(s-t53)^2+2*bi1/d^2*(s+(s-t52)^2/d)-4*bi2/d^3*(s-t52)^2+2*bi2/d^2*(s-(s-t53)^2/d);
q2=q21*(q23+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=(12*bi1/d^3-12*bi2/d^3)*s;
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q=q1+q2+q3;
end
function q=jifen_q1i_air(t51,t52,s,n,bi)
global D_d1
global L
global inter1
global eta1
d=inter1;
q11=bi/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t51))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52));
q13=1/(D_d1*(n*pi/L)^2+eta1);
q=q11*q12*q13;
end
%积分函数二，对应h4
function q=jifen_q2i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*bi1/d^2*(t52-t53)^2;
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(bi2/d^2*(t53-t52)^2);
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*bi1/d^3-12*bi2/d^3);
q44=1/(D_d1*(n*pi/L)^2+eta1);
q4=q41*(q42-q43)*q44;
q=q1+q2+q3+q4;
end
function q=jifen_q3i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*bi1/d^3-12*bi2/d^3);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q=q1+q2+q3;
end

function q=jifen_d1i_air(s,bi)
global D_d1
global L
global inter1
d=inter1;
q=bi*s;
end

function q=jifen_d2i_air(t52,t53,s,bi1,bi2)
global D_d1
global L
global inter1
d=inter1;
q1=bi1/d^2*((s-t53)^3/3+2/d*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s));
q2=bi2/d^2*((s-t52)^3/3-2/d*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s));
q=q1+q2;
end

function q = TnA_5(t,n)
global D_d1;
global eta1;
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
q111=0;q112=0;q121=0;q122=0;q131=0;q132=0;

if aa3==0
        q1=eta1*(jifen_p1i_air(t5(1),t,n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
elseif aa3==1
    if aa4==0
    q111=eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
    q112=eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
    q12=eta1*(jifen_p2i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_p2i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q13=(jifen_p3i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_p3i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q1=q111+q112+q12+q13;
    else
    q111=eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
    q112=eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
    q113=eta1*(jifen_p1i_air(t5(3),t,n,b(2))-jifen_p1i_air(t5(3),t5(3),n,b(2)));
    q121=eta1*(jifen_p2i_air(t5(2),t5(3),t5(3),n,b(1),b(2))-jifen_p2i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q122=eta1*(jifen_q2i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_q2i_air(t5(2),t5(3),t5(3),n,b(1),b(2)));
    q131=(jifen_p3i_air(t5(2),t5(3),t5(3),n,b(1),b(2))-jifen_p3i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q132=(jifen_q3i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_q3i_air(t5(2),t5(3),t5(3),n,b(1),b(2)));
    q1=q111+q112+q113+q121+q122+q131+q132;
    end
else
    if aa4==0
        q111=q111+eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
        q112=q112+eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
        for i=1:1:aa3-1
            q121=q121+eta1*(jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q122=q122+eta1*(jifen_q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q131=q131+(jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q132=q132+(jifen_q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q111=q111+eta1*(jifen_p1i_air(t5(i*2+1),t5(i*2+2),n,b(i+1))-jifen_p1i_air(t5(i*2+1),t5(i*2+1),n,b(i+1)));
            q112=q112+eta1*(jifen_q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1))-jifen_q1i_air(t5(i*2+1),t5(i*2+2),t5(i*2+2),n,b(i+1)));
        end
            q123=eta1*(jifen_p2i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1))-jifen_p2i_air(t5(aa2),t5(aa2+1),t5(aa2),n,b(aa3),b(aa3+1)));
            q133=(jifen_p3i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1))-jifen_p3i_air(t5(aa2),t5(aa2+1),t5(aa2),n,b(aa3),b(aa3+1)));
            q1=q111+q112+q121+q122+q123+q131+q132+q133;
    else
        q111=q111+eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
        q112=q112+eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
        for i=1:1:aa3-1         
            q121=q121+eta1*(jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q122=q122+eta1*(jifen_q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q131=q131+(jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q132=q132+(jifen_q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q111=q111+eta1*(jifen_p1i_air(t5(i*2+1),t5(i*2+2),n,b(i+1))-jifen_p1i_air(t5(i*2+1),t5(i*2+1),n,b(i+1)));
            q112=q112+eta1*(jifen_q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1))-jifen_q1i_air(t5(i*2+1),t5(i*2+2),t5(i*2+2),n,b(i+1)));
        end
            q121=q121+eta1*(jifen_p2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1))-jifen_p2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2),n,b(aa3),b(aa3+1)));
            q122=q122+eta1*(jifen_q2i_air(t5(aa3*2),t5(aa3*2+1),t,n,b(aa3),b(aa3+1))-jifen_q2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1)));
            q131=q131+(jifen_p3i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1))-jifen_p3i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2),n,b(aa3),b(aa3+1)));
            q132=q132+(jifen_q3i_air(t5(aa3*2),t5(aa3*2+1),t,n,b(aa3),b(aa3+1))-jifen_q3i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1)));
            q113=eta1*(jifen_p1i_air(t5(aa2),t,n,b(aa3+1))-jifen_p1i_air(t5(aa2),t5(aa2),n,b(aa3+1)));
            q1=q111+q112+q113+q121+q122+q131+q132;
    end
end






% if t>=t51 && t<t52
%     q1=eta1*(jifen_p1i_air(t51,t,n,b1)-jifen_p1i_air(t51,t51,n,b1));
% elseif t>=t52 && t<t53
%     q111=eta1*(jifen_p1i_air(t51,t52,n,b1)-jifen_p1i_air(t51,t51,n,b1));
%     q112=eta1*(jifen_q1i_air(t51,t52,t,n,b1)-jifen_q1i_air(t51,t52,t52,n,b1));
%     q12=eta1*(jifen_p2i_air(t52,t53,t,n,b1,b2)-jifen_p2i_air(t52,t53,t52,n,b1,b2));
%     q13=(jifen_p3i_air(t52,t53,t,n,b1,b2)-jifen_p3i_air(t52,t53,t52,n,b1,b2));
%     q1=q111+q112+q12+q13;
% elseif t>=t53 && t<t54
%     q111=eta1*(jifen_p1i_air(t51,t52,n,b1)-jifen_p1i_air(t51,t51,n,b1));
%     q112=eta1*(jifen_q1i_air(t51,t52,t,n,b1)-jifen_q1i_air(t51,t52,t52,n,b1));
%     q113=eta1*(jifen_p1i_air(t53,t,n,b2)-jifen_p1i_air(t53,t53,n,b2));
%     q121=eta1*(jifen_p2i_air(t52,t53,t53,n,b1,b2)-jifen_p2i_air(t52,t53,t52,n,b1,b2));
%     q122=eta1*(jifen_q2i_air(t52,t53,t,n,b1,b2)-jifen_q2i_air(t52,t53,t53,n,b1,b2));
%     q131=(jifen_p3i_air(t52,t53,t53,n,b1,b2)-jifen_p3i_air(t52,t53,t52,n,b1,b2));
%     q132=(jifen_q3i_air(t52,t53,t,n,b1,b2)-jifen_q3i_air(t52,t53,t53,n,b1,b2));
%     q1=q111+q112+q113+q121+q122+q131+q132;
% elseif t>=t54 && t<t55
%     q111=eta1*(jifen_p1i_air(t51,t52,n,b1)-jifen_p1i_air(t51,t51,n,b1));
%     q112=eta1*(jifen_q1i_air(t51,t52,t,n,b1)-jifen_q1i_air(t51,t52,t52,n,b1));
%     q113=eta1*(jifen_p1i_air(t53,t54,n,b2)-jifen_p1i_air(t53,t53,n,b2));
%     q114=eta1*(jifen_q1i_air(t53,t54,t,n,b2)-jifen_q1i_air(t53,t54,t54,n,b2));
%     q121=eta1*(jifen_p2i_air(t52,t53,t53,n,b1,b2)-jifen_p2i_air(t52,t53,t52,n,b1,b2));
%     q122=eta1*(jifen_q2i_air(t52,t53,t,n,b1,b2)-jifen_q2i_air(t52,t53,t53,n,b1,b2));
%     q123=eta1*(jifen_p2i_air(t54,t55,t,n,b2,b3)-jifen_p2i_air(t54,t55,t54,n,b2,b3));
%     q131=(jifen_p3i_air(t52,t53,t53,n,b1,b2)-jifen_p3i_air(t52,t53,t52,n,b1,b2));
%     q132=(jifen_q3i_air(t52,t53,t,n,b1,b2)-jifen_q3i_air(t52,t53,t53,n,b1,b2));
%     q133=(jifen_p3i_air(t54,t55,t,n,b2,b3)-jifen_p3i_air(t54,t55,t54,n,b2,b3));
%     q1=q111+q112+q113+q114+q121+q122+q123+q131+q132+q133;
% elseif t>=t55 && t<t56
%     q111=eta1*(jifen_p1i_air(t51,t52,n,b1)-jifen_p1i_air(t51,t51,n,b1));
%     q112=eta1*(jifen_q1i_air(t51,t52,t,n,b1)-jifen_q1i_air(t51,t52,t52,n,b1));
%     q113=eta1*(jifen_p1i_air(t53,t54,n,b2)-jifen_p1i_air(t53,t53,n,b2));
%     q114=eta1*(jifen_q1i_air(t53,t54,t,n,b2)-jifen_q1i_air(t53,t54,t54,n,b2));
%     q115=eta1*(jifen_p1i_air(t55,t,n,b3)-jifen_p1i_air(t55,t55,n,b3));
%     q121=eta1*(jifen_p2i_air(t52,t53,t53,n,b1,b2)-jifen_p2i_air(t52,t53,t52,n,b1,b2));
%     q122=eta1*(jifen_q2i_air(t52,t53,t,n,b1,b2)-jifen_q2i_air(t52,t53,t53,n,b1,b2));
%     q123=eta1*(jifen_p2i_air(t54,t55,t55,n,b2,b3)-jifen_p2i_air(t54,t55,t54,n,b2,b3));
%     q124=eta1*(jifen_q2i_air(t54,t55,t,n,b2,b3)-jifen_q2i_air(t54,t55,t55,n,b2,b3));
%     q131=(jifen_p3i_air(t52,t53,t53,n,b1,b2)-jifen_p3i_air(t52,t53,t52,n,b1,b2));
%     q132=(jifen_q3i_air(t52,t53,t,n,b1,b2)-jifen_q3i_air(t52,t53,t53,n,b1,b2));
%     q133=(jifen_p3i_air(t54,t55,t55,n,b2,b3)-jifen_p3i_air(t54,t55,t54,n,b2,b3));
%     q134=(jifen_q3i_air(t54,t55,t,n,b2,b3)-jifen_q3i_air(t54,t55,t55,n,b2,b3));
%     q1=q111+q112+q113+q114+q115+q121+q122+q123+q124+q131+q132+q133+q134;
% elseif t>=t56 && t<t57
%     q111=eta1*(jifen_p1i_air(t51,t52,n,b1)-jifen_p1i_air(t51,t51,n,b1));
%     q112=eta1*(jifen_q1i_air(t51,t52,t,n,b1)-jifen_q1i_air(t51,t52,t52,n,b1));
%     q113=eta1*(jifen_p1i_air(t53,t54,n,b2)-jifen_p1i_air(t53,t53,n,b2));
%     q114=eta1*(jifen_q1i_air(t53,t54,t,n,b2)-jifen_q1i_air(t53,t54,t54,n,b2));
%     q115=eta1*(jifen_p1i_air(t55,t56,n,b3)-jifen_p1i_air(t55,t55,n,b3));
%     q116=eta1*(jifen_q1i_air(t55,t56,t,n,b3)-jifen_q1i_air(t55,t56,t56,n,b3));
%     q121=eta1*(jifen_p2i_air(t52,t53,t53,n,b1,b2)-jifen_p2i_air(t52,t53,t52,n,b1,b2));
%     q122=eta1*(jifen_q2i_air(t52,t53,t,n,b1,b2)-jifen_q2i_air(t52,t53,t53,n,b1,b2));
%     q123=eta1*(jifen_p2i_air(t54,t55,t55,n,b2,b3)-jifen_p2i_air(t54,t55,t54,n,b2,b3));
%     q124=eta1*(jifen_q2i_air(t54,t55,t,n,b2,b3)-jifen_q2i_air(t54,t55,t55,n,b2,b3));
%     q125=eta1*(jifen_p2i_air(t56,t57,t,n,b3,b4)-jifen_p2i_air(t56,t57,t56,n,b3,b4));
%     q131=(jifen_p3i_air(t52,t53,t53,n,b1,b2)-jifen_p3i_air(t52,t53,t52,n,b1,b2));
%     q132=(jifen_q3i_air(t52,t53,t,n,b1,b2)-jifen_q3i_air(t52,t53,t53,n,b1,b2));
%     q133=(jifen_p3i_air(t54,t55,t55,n,b2,b3)-jifen_p3i_air(t54,t55,t54,n,b2,b3));
%     q134=(jifen_q3i_air(t54,t55,t,n,b2,b3)-jifen_q3i_air(t54,t55,t55,n,b2,b3));
%     q135=(jifen_p3i_air(t56,t57,t,n,b3,b4)-jifen_p3i_air(t56,t57,t56,n,b3,b4));
%     q1=q111+q112+q113+q114+q115+q116+q121+q122+q123+q124+q125+q131+q132+q133+q134+q135;
% elseif t>=t57 && t<=t58
%     q111=eta1*(jifen_p1i_air(t51,t52,n,b1)-jifen_p1i_air(t51,t51,n,b1));
%     q112=eta1*(jifen_q1i_air(t51,t52,t,n,b1)-jifen_q1i_air(t51,t52,t52,n,b1));
%     q113=eta1*(jifen_p1i_air(t53,t54,n,b2)-jifen_p1i_air(t53,t53,n,b2));
%     q114=eta1*(jifen_q1i_air(t53,t54,t,n,b2)-jifen_q1i_air(t53,t54,t54,n,b2));
%     q115=eta1*(jifen_p1i_air(t55,t56,n,b3)-jifen_p1i_air(t55,t55,n,b3));
%     q116=eta1*(jifen_q1i_air(t55,t56,t,n,b3)-jifen_q1i_air(t55,t56,t56,n,b3));
%     q117=eta1*(jifen_p1i_air(t57,t,n,b4)-jifen_p1i_air(t57,t57,n,b4));
%     q121=eta1*(jifen_p2i_air(t52,t53,t53,n,b1,b2)-jifen_p2i_air(t52,t53,t52,n,b1,b2));
%     q122=eta1*(jifen_q2i_air(t52,t53,t,n,b1,b2)-jifen_q2i_air(t52,t53,t53,n,b1,b2));
%     q123=eta1*(jifen_p2i_air(t54,t55,t55,n,b2,b3)-jifen_p2i_air(t54,t55,t54,n,b2,b3));
%     q124=eta1*(jifen_q2i_air(t54,t55,t,n,b2,b3)-jifen_q2i_air(t54,t55,t55,n,b2,b3));
%     q125=eta1*(jifen_p2i_air(t56,t57,t57,n,b3,b4)-jifen_p2i_air(t56,t57,t56,n,b3,b4));
%     q126=eta1*(jifen_q2i_air(t56,t57,t,n,b3,b4)-jifen_q2i_air(t56,t57,t57,n,b3,b4));
%     q131=(jifen_p3i_air(t52,t53,t53,n,b1,b2)-jifen_p3i_air(t52,t53,t52,n,b1,b2));
%     q132=(jifen_q3i_air(t52,t53,t,n,b1,b2)-jifen_q3i_air(t52,t53,t53,n,b1,b2));
%     q133=(jifen_p3i_air(t54,t55,t55,n,b2,b3)-jifen_p3i_air(t54,t55,t54,n,b2,b3));
%     q134=(jifen_q3i_air(t54,t55,t,n,b2,b3)-jifen_q3i_air(t54,t55,t55,n,b2,b3));
%     q135=(jifen_p3i_air(t56,t57,t57,n,b3,b4)-jifen_p3i_air(t56,t57,t56,n,b3,b4));
%     q136=(jifen_q3i_air(t56,t57,t,n,b3,b4)-jifen_q3i_air(t56,t57,t57,n,b3,b4));
%     q1=q111+q112+q113+q114+q115+q116+q117+q121+q122+q123+q124+q125+q126+q131+q132+q133+q134+q135+q136;
% else
%     msg='Error occured';
%     error(msg)
% end

fi_n=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n);
q4=2*L/(n*pi)*(-1)^n*q1;
q=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1)+q4;
end

%%这个是用药为阶梯l两段函数，一段时间用药b1，一段时间用药b2，一段时间用药b3,一段时间用药b4，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q = A_5(x,t) 
global L
global t51 t52 t53 t54 t55 t56 t57 t58 
global b1 b2 b3 b4
global A_0
global B_3  %A的分泌速率 
global b t5

aa1=find(t>=t5);
aa2=aa1(end);
aa3=floor(aa2/2);
aa4=aa2-aa3*2;

q11=0;q12=0;
if aa3==0
        q1=jifen_d1i_air(t,b(1))-jifen_d1i_air(t5(1),b(1));
elseif aa3==1
    if aa4==0
        q11=jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        q12=jifen_d2i_air(t5(2),t5(3),t,b(1),b(2))-jifen_d2i_air(t5(2),t5(3),t5(2),b(1),b(2));
        q1=q11+q12;
    else
        q11=jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        q12=jifen_d2i_air(t5(2),t5(3),t5(3),b(1),b(2))-jifen_d2i_air(t5(2),t5(3),t5(2),b(1),b(2));
        q13=jifen_d1i_air(t,b(2))-jifen_d1i_air(t5(3),b(2));
        q1=q11+q12+q13;
    end
else
    if aa4==0
        q11=q11+jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        for i=1:1:aa3-1
            q12=q12+jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),b(i),b(i+1))-jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2),b(i),b(i+1));
            q11=q11+jifen_d1i_air(t5(i*2+2),b(i+1))-jifen_d1i_air(t5(i*2+1),b(i+1));
        end
        q12=q12+jifen_d2i_air(t5(aa2),t5(aa2+1),t,b(aa3),b(aa3+1))-jifen_d2i_air(t5(aa2),t5(aa2+1),t5(aa2),b(aa3),b(aa3+1));
        q1=q11+q12;
    else
        q11=q11+jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        for i=1:1:aa3-1
            q12=q12+jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),b(i),b(i+1))-jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2),b(i),b(i+1));
            q11=q11+jifen_d1i_air(t5(i*2+2),b(i+1))-jifen_d1i_air(t5(i*2+1),b(i+1));
        end
        q12=q12+jifen_d2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),b(aa3),b(aa3+1))-jifen_d2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2),b(aa3),b(aa3+1));
        q11=q11+jifen_d1i_air(t,b(aa3+1))-jifen_d1i_air(t5(aa2),b(aa3+1));
        q1=q11+q12;
    end
end
dnt=q1;



% if t>=t51 && t<t52
%     dnt=jifen_d1i_air(t,b1)-jifen_d1i_air(t51,b1);
% elseif t>=t52 && t<t53
%     dnt1=jifen_d1i_air(t52,b1)-jifen_d1i_air(t51,b1);
%     dnt2=jifen_d2i_air(t52,t53,t,b1,b2)-jifen_d2i_air(t52,t53,t52,b1,b2);
%     dnt=dnt1+dnt2;
% elseif t>=t53 && t<t54
%     dnt1=jifen_d1i_air(t52,b1)-jifen_d1i_air(t51,b1);
%     dnt2=jifen_d2i_air(t52,t53,t53,b1,b2)-jifen_d2i_air(t52,t53,t52,b1,b2);
%     dnt3=jifen_d1i_air(t,b2)-jifen_d1i_air(t53,b2);
%     dnt=dnt1+dnt2+dnt3;
% elseif t>=t54 && t<t55
%     dnt1=jifen_d1i_air(t52,b1)-jifen_d1i_air(t51,b1);
%     dnt2=jifen_d2i_air(t52,t53,t53,b1,b2)-jifen_d2i_air(t52,t53,t52,b1,b2);
%     dnt3=jifen_d1i_air(t54,b2)-jifen_d1i_air(t53,b2);
%     dnt4=jifen_d2i_air(t54,t55,t,b2,b3)-jifen_d2i_air(t54,t55,t54,b2,b3);
%     dnt=dnt1+dnt2+dnt3+dnt4;
% elseif t>=t55 && t<t56
%     dnt1=jifen_d1i_air(t52,b1)-jifen_d1i_air(t51,b1);
%     dnt2=jifen_d2i_air(t52,t53,t53,b1,b2)-jifen_d2i_air(t52,t53,t52,b1,b2);
%     dnt3=jifen_d1i_air(t54,b2)-jifen_d1i_air(t53,b2);
%     dnt4=jifen_d2i_air(t54,t55,t55,b2,b3)-jifen_d2i_air(t54,t55,t54,b2,b3);
%     dnt5=jifen_d1i_air(t,b3)-jifen_d1i_air(t55,b3);
%     dnt=dnt1+dnt2+dnt3+dnt4+dnt5;
% elseif t>=t56 && t<t57
%     dnt1=jifen_d1i_air(t52,b1)-jifen_d1i_air(t51,b1);
%     dnt2=jifen_d2i_air(t52,t53,t53,b1,b2)-jifen_d2i_air(t52,t53,t52,b1,b2);
%     dnt3=jifen_d1i_air(t54,b2)-jifen_d1i_air(t53,b2);
%     dnt4=jifen_d2i_air(t54,t55,t55,b2,b3)-jifen_d2i_air(t54,t55,t54,b2,b3);
%     dnt5=jifen_d1i_air(t56,b3)-jifen_d1i_air(t55,b3);
%     dnt6=jifen_d2i_air(t56,t57,t,b3,b4)-jifen_d2i_air(t56,t57,t56,b3,b4);
%     dnt=dnt1+dnt2+dnt3+dnt4+dnt5+dnt6;
% elseif t>=t57 && t<=t58
%     dnt1=jifen_d1i_air(t52,b1)-jifen_d1i_air(t51,b1);
%     dnt2=jifen_d2i_air(t52,t53,t53,b1,b2)-jifen_d2i_air(t52,t53,t52,b1,b2);
%     dnt3=jifen_d1i_air(t54,b2)-jifen_d1i_air(t53,b2);
%     dnt4=jifen_d2i_air(t54,t55,t55,b2,b3)-jifen_d2i_air(t54,t55,t54,b2,b3);
%     dnt5=jifen_d1i_air(t56,b3)-jifen_d1i_air(t55,b3);
%     dnt6=jifen_d2i_air(t56,t57,t57,b3,b4)-jifen_d2i_air(t56,t57,t56,b3,b4);
%     dnt7=jifen_d1i_air(t,b4)-jifen_d1i_air(t57,b4);
%     dnt=dnt1+dnt2+dnt3+dnt4+dnt5+dnt6+dnt7;
% % if t>=t51 && t<t52
% %     dnt=b1;
% % elseif t>=t52 && t<t53
% %     dnt=b1/d^2*(1+2/d*(t-t52))*(t-t53)^2+b2/d^2*(1-2/d*(t-t53))*(t-t52)^2;
% % elseif t>=t53
% %     dnt=b2;
% % else
% else
%     msg='Error occured';
%     error(msg)
% end

q0=0;
    for n=1:1:20
    q0=q0+TnA_5(t,n)*sin(n*pi*x/L);  
    end
q=A_0+B_3/x*q0+B_3*dnt;
end