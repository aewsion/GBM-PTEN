function aaa = Drug(x,t,select,dmax,dmax_1)
% global eta1    %药物的降解率
% global D_d1    %药物的扩散率
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
        [aaa bbb]=picture(x,t);
    case 2
        [bbb aaa]=picture(x,t);
    case 3
        aaa=Drug_34(x,t);
    case 4
        aaa=Drug_ceshi_jietihanshu(x,t);
    case 5
        aaa=Drug_5(x,t) ;
end
end
function [drugt1,drugt2] = picture(x,t)
global d_max
global eta1
global D_d1
global L   %x的最大值
global delt_t
global TimeLength period1
d_initial=d_max/L;
w=pi/period1;
% w=pi/(TimeLength*delt_t/(period1/delt_t));
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

%%这个是用药为阶梯函数，一段时间用药6，一段时间用药0，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第一部分
%积分函数一，对应h1(t,t1)  t1是to
function q4=h_1_ceshiair1(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*exp(-D_d1*(n*pi/L)^2*(t-t1))*(2*m0/d^3*(t1-(k+1/2)*T)^2+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d))*(t1-(k+1/2)*T));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*exp(-D_d1*(n*pi/L)^2*(t-t1))*(8*m0/d^3*(t1-(k+1/2)*T)+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d)));
q3=24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(t-t1));
q4=q1+q2+q3;
end
%积分函数二，对应h1(t,t1)
function q4=h_2_ceshiair2(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=2*(-1)^n/D_d1/(n*pi/L)^3*exp(-D_d1*(n*pi/L)^2*(t-t1))*(-2*m0/d^3*(t1-(k+1)*T+d)^2+2*m0/d^2*(1-2/d*(t1-(k+1)*T))*(t1-(k+1)*T+d));
q2=-2*(-1)^n/D_d1^2/(n*pi/L)^5*exp(-D_d1*(n*pi/L)^2*(t-t1))*(-8*m0/d^3*(t1-(k+1)*T+d)+2*m0/d^2*(1-2/d*(t1-(k+1)*T)));
q3=-24*m0*(-1)^n/d^3/D_d1^3/(n*pi/L)^7*exp(-D_d1*(n*pi/L)^2*(t-t1));
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
global m0
q0=0;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
xk=(d1+1)*period1-inter1;
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

%下面是测试函数阶梯函数
function Drug4 = Drug_ceshi_jietihanshu(x,t)
global d_max
global eta1
global D_d1
global L   %x的最大值
global delt_t
global TimeLength period1
d_initial=d_max/L;
w=pi/period1;
% w=pi/(TimeLength*delt_t/(period1/delt_t));
r1=floor(t/period1);
r2=rem(r1,2);
if r2==0
    d11=d_max;
else
    d11=0;
end
%    d11=d_max;
   d110=d_max;
drug4=0;
for n=1:1:100
    Tn1=((2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+  d110*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t)))*sin(n*pi*x/L);
    drug4=drug4+Tn1; 
end
Drug4=1/x*drug4+d11;
end

%下面是用药停药函数的函数，药物情况5
function q=h_3_air(t52,t53,t,n)
global D_d1
global L
global inter1
global period1
global b1 b2
d=inter1;
q1=2*b1/d^3*L^2/(D_d1*(n*pi)^2)*((t-t53)^2-exp(D_d1*(n*pi/L)^2*(t52-t))*(t52-t53)^2);
q2=-4*b1/d^3*L^4/(D_d1^2*(n*pi)^4)*((t-t53)-exp(D_d1*(n*pi/L)^2*(t52-t))*(t52-t53));
q3=4*b1/d^3*L^6/(D_d1^3*(n*pi)^6)*(1-exp(D_d1*(n*pi/L)^2*(t52-t)));
q4=2*b1/d^2*L^2/(D_d1*(n*pi)^2)*((1+2/d*(t-t52))*(t-t53)-exp(D_d1*(n*pi/L)^2*(t52-t))*(t52-t53));
q5=-2*b1/d^2*L^4/(D_d1^2*(n*pi)^4)*(1+2/d*(2*t-t53-t52)-exp(D_d1*(n*pi/L)^2*(t52-t))*(2/d*(t52-t53)+1));
q6=8*b1/d^3*L^6/(D_d1^3*(n*pi)^6)*(1-exp(D_d1*(n*pi/L)^2*(t52-t)));
q7=-2*b2/d^3*L^2/(D_d1*(n*pi)^2)*(t-t52)^2;
q8=4*b2/d^3*L^4/(D_d1^2*(n*pi)^4)*(t-t52);
q9=-4*b2/d^3*L^6/(D_d1^3*(n*pi)^6)*(1-exp(D_d1*(n*pi/L)^2*(t52-t)));
q10=2*b2/d^2*L^2/(D_d1*(n*pi)^2)*(1-2/d*(t-t53))*(t-t52);
q11=-2*b2/d^2*L^4/(D_d1^2*(n*pi)^4)*(1-2/d*(2*t-t52-t53)-exp(D_d1*(n*pi/L)^2*(t52-t))*(1-2/d*(t52-t53)));
q12=-8*b2/d^3*L^6/(D_d1^3*(n*pi)^6)*(1-exp(D_d1*(n*pi/L)^2*(t52-t)));
q=q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12;
end
%积分函数二，对应h4
function q=h_4_air(t52,t53,t,n)
global D_d1
global L
global inter1
global period1
global b1 b2
d=inter1;
q1=2*b1/d^3*L^2/(D_d1*(n*pi)^2)*(-exp(D_d1*(n*pi/L)^2*(t52-t))*(t52-t53)^2);
q2=-4*b1/d^3*L^4/(D_d1^2*(n*pi)^4)*(-exp(D_d1*(n*pi/L)^2*(t52-t))*(t52-t53));
q3=4*b1/d^3*L^6/(D_d1^3*(n*pi)^6)*(exp(D_d1*(n*pi/L)^2*(t53-t))-exp(D_d1*(n*pi/L)^2*(t52-t)));
q4=2*b1/d^2*L^2/(D_d1*(n*pi)^2)*(-exp(D_d1*(n*pi/L)^2*(t52-t))*(t52-t53));
q5=-2*b1/d^2*L^4/(D_d1^2*(n*pi)^4)*(exp(D_d1*(n*pi/L)^2*(t53-t))*(1+2/d*(t53-t52))-exp(D_d1*(n*pi/L)^2*(t52-t))*(2/d*(t52-t53)+1));
q6=8*b1/d^3*L^6/(D_d1^3*(n*pi)^6)*(exp(D_d1*(n*pi/L)^2*(t53-t))-exp(D_d1*(n*pi/L)^2*(t52-t)));
q7=-2*b2/d^3*L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t53-t))*(t53-t52)^2;
q8=4*b2/d^3*L^4/(D_d1^2*(n*pi)^4)*exp(D_d1*(n*pi/L)^2*(t53-t))*(t53-t52);
q9=-4*b2/d^3*L^6/(D_d1^3*(n*pi)^6)*(exp(D_d1*(n*pi/L)^2*(t53-t))-exp(D_d1*(n*pi/L)^2*(t52-t)));
q10=2*b2/d^2*L^2/(D_d1*(n*pi)^2)*exp(D_d1*(n*pi/L)^2*(t53-t))*(t53-t52);
q11=-2*b2/d^2*L^4/(D_d1^2*(n*pi)^4)*(exp(D_d1*(n*pi/L)^2*(t53-t))*(1-2/d*(t53-t52))-exp(D_d1*(n*pi/L)^2*(t52-t))*(1-2/d*(t52-t53)));
q12=-8*b2/d^3*L^6/(D_d1^3*(n*pi)^6)*(exp(D_d1*(n*pi/L)^2*(t53-t))-exp(D_d1*(n*pi/L)^2*(t52-t)));
q=q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12;
end
%%这个是用药为阶梯函数，一段时间用药b1，一段时间用药b2，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第二部分
function q = Tn5(t,n)
global D_d1;
global eta1;
global L
global period1
global inter1
global b1 b2
global t51 t52 t53
d_0=b1;%初始条件为d_initial*r，应该要除以L
d_initial=d_0/L;
if t>=t51 && t<t52
    q1=0;
elseif t>=t52 && t<t53
    q1=h_3_air(t52,t53,t,n);
elseif t>=t53
    q1=h_4_air(t52,t53,t,n);
else
    msg='Error occured';
    error(msg)
end
q2=q1*2*L/(n*pi)*(-1)^n;
q=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n)*exp(-D_d1*(n*pi/L)^2*t)+eta1*2*L^3/(n*pi)^3/D_d1*(-1)^n*(1-exp(-D_d1*(n*pi/L)^2*t))+q2;
end


%%这个是用药为阶梯l两段函数，一段时间用药b1，一段时间用药b2，连接函数用埃尔米插值来连接，用来计算Tn,其中涉及到积分就用累加的积分
%第三部分
function q = Drug_5(x,t) 
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
    q0=q0+Tn5(t,n)*sin(n*pi*x/L);  
    end
q=q0/x+dnt;
end