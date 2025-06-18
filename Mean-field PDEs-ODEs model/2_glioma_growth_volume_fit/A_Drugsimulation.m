function aaa=A_Drugsimulation(x,t,select,dmax,dmax_1)
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
    b1=dmax;
% end
if nargin >=5
    b2=dmax_1;
end
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
function q = TnA_1(x,t,n)
global d_max
global D_d1;
global eta1;
global L
d_initial=d_max/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
% q2=0;%因为又边界输入为常数3，导数值为0；
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_max*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/(D_d1*(n*pi/L)^2)+eta1*2*(-1)^n*t/(D_d1*(n*pi/L)^3)+eta1*2*(-1)^n*(exp(-D_d1*(n*pi/L)^2*t)-1)/(D_d1^2*(n*pi/L)^5);

end
%第二部分
function q = A_1(x,t)
global d_max
global A_0
global L
global B_3  %A的分泌速率 
q0=0;

for n=1:1:100
    q0=q0+sin(n*pi/L*x)*TnA_1(x,t,n);
end
q=A_0+B_3/x*q0+B_3*d_max*t;
end

%%药物边界条件为d22=d_max/2*(sin(w*t)+2)=3/2*(sin(w*t)+2)的时候产生的A的量
%第一部分
function q5 = TnA_2(t,n)
global d_max
global D_d1;
global eta1;
global L
global delt_t
global TimeLength period1
% w=pi/(TimeLength*delt_t/20);
w=pi/period1;
d_initial=d_max/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
% q2=0;%因为又边界输入为常数3，导数值为0；
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_max*(2*L)/n/pi*(-1)^n;
q1=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/(D_d1*(n*pi/L)^2);
q2=eta1*2*(-1)^n*t/(D_d1*(n*pi/L)^3);
q3=eta1*2*(-1)^n/(D_d1^2*(n*pi/L)^5)*(exp(-D_d1*(n*pi/L)^2*t)-1);
q4=d_max*w*L/n/pi*(-1)^n*(D_d1*(n*pi/L)^2*sin(w*t)/w+exp(-D_d1*(n*pi/L)^2*t)-cos(w*t))/(w^2+D_d1^2*(n*pi/L)^4);
q5=q1+q2+q3+q4;
end
%第二部分
function q = A_2(x,t)
global d_max
global A_0
global L
global B_3  %A的分泌速率 
global TimeLength
global delt_t
q0=0;
w=pi/(TimeLength*delt_t/20);
for n=1:1:100
    q0=q0+sin(n*pi*x/L)*TnA_2(t,n);
end
q=A_0+B_3/x*q0+B_3*d_max/2/w*(1-cos(w*t))+B_3*d_max*t;
end

%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分


%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数三
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
%第一部分
function q6=h_3_ceshiair3(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*(2*m0/3/d^3*(s-(k+1/2)*T)^3+2*m0/d^2*(s^2/2+2/d*(s^3/3-((k+1/2)*T-d)/2*s^2))-2*m0*(k+1/2)*T/d^2*(s+1/d*(s-(k+1/2)*T+d)^2));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*(4*m0/d^3*(s-(k+1/2)*T)^2+2*m0/d^2*(s+1/d*(s-(k+1/2)*T+d)^2));
q3=24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*s;
q4=12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q5=24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q6=q1+q2+q3+q4+q5;
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数四
function q6=h_4_ceshiair4(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*(-2*m0/3/d^3*(s-(k+1)*T+d)^3+2*m0/d^2*(s^2/2-2/d*(s^3/3-(k+1)*T/2*s^2))-2*m0*((k+1)*T-d)/d^2*(s-1/d*(s-(k+1)*T)^2));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*(-4*m0/d^3*(s-(k+1)*T+d)^2+2*m0/d^2*(s-1/d*(s-(k+1)*T)^2));
q3=-24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*s;
q4=-12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q5=-24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q6=q1+q2+q3+q4+q5;
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数五
%对应的是对gk1函数的积分
%积分函数五是因为fn包含边界函数d'(t)是阶梯函数，分段积分完成了后是一个定的关于t的积分，再对t积分后和积分函数三不一样，然后才是IL的值的一部分
function q5=h_5_ceshiair5(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T));
q2=12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q3=-24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T));
q4=24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1/2)*T+d));
q5=q1+q2+q3+q4;
end
%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为积分函数六
%对应的是对gk2函数的积分
function q5=h_6_ceshiair6(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=-12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T));
q2=-12*m0*(-1)^n/d^2/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q3=24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T));
q4=-24*m0*(-1)^n/d^3/D_d1^4/(n*pi/L)^9*exp(-D_d1*(n*pi/L)^2*(s-(k+1)*T+d));
q5=q1+q2+q3+q4;
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
global D_d1
global L
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
global D_d1
global L
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
    q0=q0+TnA_3(t,n)*sin(n*pi*x/L);
end
% q=A_0+B_3/x*q0+B_3*d_0*t+B_3*q4;
q=A_0+B_3/x*q0+B_3*q4;
end

function q = TnA_4(x,t,n)
global d_max
global D_d1;
global eta1;
global L
d_initial=d_max/L;%初始条件为d_initial*r--x的意思，所以应该要除以100
% q2=0;%因为又边界输入为常数3，导数值为0；
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_max*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/(D_d1*(n*pi/L)^2)+eta1*2*(-1)^n*t/(D_d1*(n*pi/L)^3)+eta1*2*(-1)^n*(exp(-D_d1*(n*pi/L)^2*t)-1)/(D_d1^2*(n*pi/L)^5);

end

%第二部分
function q = A_4(x,t)
global d_max
global A_0
global L
global B_3  %A的分泌速率 
global period1  %period1为一个周期

r1=floor(t/period1);
r2=rem(r1,2);
r3=floor(r1/2);
if r2==0
    d21=d_max*(t-r1*period1);
else
    d21=d_max*period1;
end
d22=r3*period1*d_max;

q0=0;

for n=1:1:100
    q0=q0+sin(n*pi/L*x)*TnA_4(x,t,n);
end
% q=A_0+B_3/x*q0+B_3*d_max*t;
%新的阶梯函数的
q=A_0+B_3/x*q0+B_3*(d21+d22);
end

%%用药情况5：这个是用药为阶梯分两段函数，一段时间用药b1，一段时间用药b2，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分

%用来对于阶梯函数用埃尔米插值法的时候算A，对药物的一部分进行积分，下面为函数三的积分函数
%积分函数三是因为fn包含边界函数d'(t)是阶梯函数，要分段积分，积分后再对t积分才是IL的值的一部分
function q=h_3_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global b1 b2
d=inter1;
q1=2*b1/d^3*L^2/(D_d1*(n*pi)^2)*((s-t53)^3/3+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s))*(t52-t53)^2);
q2=-4*b1/d^3*L^4/(D_d1^2*(n*pi)^4)*((s-t53)^2/2+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s))*(t52-t53));
q3=4*b1/d^3*L^6/(D_d1^3*(n*pi)^6)*(s+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s)));
q4=2*b1/d^2*L^2/(D_d1*(n*pi)^2)*((s-t53)^2/2+2/d*(s^3/3-t52*s^2/2-t53*s^2/2+t52*t53*s)+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s))*(t52-t53));
q5=-2*b1/d^2*L^4/(D_d1^2*(n*pi)^4)*(s+2/d*(s^2-t53*s-t52*s)+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s))*(2/d*(t52-t53)+1));
q6=8*b1/d^3*L^6/(D_d1^3*(n*pi)^6)*(s+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s)));
q7=-2*b2/d^3*L^2/(D_d1*(n*pi)^2)*(s-t52)^3/3;
q8=4*b2/d^3*L^4/(D_d1^2*(n*pi)^4)*(s-t52)^2/2;
q9=-4*b2/d^3*L^6/(D_d1^3*(n*pi)^6)*(s+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s)));
q10=2*b2/d^2*L^2/(D_d1*(n*pi)^2)*((s-t52)^2/2-2/d*(s^3/3-t53*s^2/2-t52*s^2/2+t52*t53*s));
q11=-2*b2/d^2*L^4/(D_d1^2*(n*pi)^4)*(s-2/d*(s^2-t52*s-t53*s)+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s))*(1-2/d*(t52-t53)));
q12=-8*b2/d^3*L^6/(D_d1^3*(n*pi)^6)*(s+L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t52-s)));
q=q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12;
end
function q=h_4_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global b1 b2
d=inter1;
q1=2*b1/d^3*L^4/(D_d1^2*(n*pi)^4)*exp(D_d1*(n*pi/L)^2*(t52-s))*(t52-t53)^2;
q2=-4*b1/d^3*L^6/(D_d1^3*(n*pi)^6)*exp(D_d1*(n*pi/L)^2*(t52-s))*(t52-t53);
q3=4*b1/d^3*L^8/(D_d1^4*(n*pi)^8)*(exp(D_d1*(n*pi/L)^2*(t52-s))-exp(D_d1*(n*pi/L)^2*(t53-s)));
q4=2*b1/d^2*L^4/(D_d1^2*(n*pi)^4)*exp(D_d1*(n*pi/L)^2*(t52-s))*(t52-t53);
q5=-2*b1/d^2*L^6/(D_d1^3*(n*pi)^6)*(exp(D_d1*(n*pi/L)^2*(t52-s))*(2/d*(t52-t53)+1)-exp(D_d1*(n*pi/L)^2*(t53-s))*(1+2/d*(t53-t52)));
q6=8*b1/d^3*L^8/(D_d1^4*(n*pi)^8)*(exp(D_d1*(n*pi/L)^2*(t52-s))-exp(D_d1*(n*pi/L)^2*(t53-s)));
q7=2*b2/d^3*L^4/(D_d1^2*(n*pi)^4)*exp(D_d1*(n*pi/L)^2*(t53-s))*(t53-t52)^2;
q8=-4*b2/d^3*L^6/(D_d1^3*(n*pi)^6)*exp(D_d1*(n*pi/L)^2*(t53-s))*(t53-t52);
q9=-4*b2/d^3*L^8/(D_d1^4*(n*pi)^8)*(exp(D_d1*(n*pi/L)^2*(t52-s))-exp(D_d1*(n*pi/L)^2*(t53-s)));
q10=-2*b2/d^2*L^4/(D_d1^2*(n*pi)^4)*exp(D_d1*(n*pi/L)^2*(t53-s))*(t53-t52);
q11=-2*b2/d^2*L^6/(D_d1^3*(n*pi)^6)*(exp(D_d1*(n*pi/L)^2*(t52-s))*(1-2/d*(t52-t53))-exp(D_d1*(n*pi/L)^2*(t53-s))*(1-2/d*(t53-t52)));
q12=-8*b2/d^3*L^8/(D_d1^4*(n*pi)^8)*(exp(D_d1*(n*pi/L)^2*(t52-s))-exp(D_d1*(n*pi/L)^2*(t53-s)));
q=q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12;
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
function q = TnA_5(t,n)
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
    q4=0;
elseif t>=t52 && t<t53
    q41=h_3_airjifen(t52,t53,t52,n);
    q42=h_3_airjifen(t52,t53,t,n);
    q4=q42-q41;
elseif t>=t53
    q41=h_3_airjifen(t52,t53,t52,n);
    q42=h_3_airjifen(t52,t53,t53,n);
    q43=h_4_airjifen(t52,t53,t53,n);
    q44=h_4_airjifen(t52,t53,t,n);
    q4=q42-q41+q44-q43;
else
    msg='Error occured';
    error(msg)    
end
q5=q4*2*L/n/pi*(-1)^n;
fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q1=fi_n*(1-exp(-D_d1*(n*pi/L)^2*t))/D_d1/(n*pi/L)^2;
q2=eta1*2*(-1)^n/D_d1/(n*pi/L)^3*t;
q3=eta1*2*(-1)^n/D_d1^2/(n*pi/L)^5*(exp(-D_d1*(n*pi/L)^2*t)-1);
q=q1+q2+q3+q5;
end
%第三部分
function q = A_5(x,t)
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
    q0=q0+TnA_5(t,n)*sin(n*pi*x/L);
end
% q=A_0+B_3/x*q0+B_3*d_0*t+B_3*q4;
q=A_0+B_3/x*q0+B_3*q4;
end