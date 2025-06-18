function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
global period1 d_max period_2 delt_t L 

% d11=1;
% 
% C_d=d_max;
% w=2*pi/period_2;
% d22=C_d/2*(sin(w*t)+2);




inter1=1*delt_t
m0=1;
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
d33=dnt;

% global L
% global t51 t52 t53
% global inter1
% global b1 b2
% d=inter1;
% if t>=t51 && t<t52
%     dnt=b1;
% elseif t>=t52 && t<t53
%     dnt=b1/d^2*(1+2/d*(t-t52))*(t-t53)^2+b2/d^2*(1-2/d*(t-t53))*(t-t52)^2;
% elseif t>=t53
%     dnt=b2;
% else
%         msg='Error occured';
%     error(msg)
% end
% d44=dnt;

pl=[0];
ql=[1];
% pr=[ur-d11];
% pr=[ur-d22];
pr=[ur-d33];
% pr=[ur-d44];
% pr=[ur-d55];
qr=[0];

%pr = [0; 0; 0; 0; ur(5)-dnt];
%pr = [0; 0; 0; 0; ur(5)-d_max/2*(sin((t/delt_t)*pi/100)+2)];
%pr = [0; 0; 0; 0; ur(5)-(3-rem(floor(t/(delt_t*30)),3))];%忘记了是哪种情况了
% qr = [1; 1; 1; 1; 1; 1];
% pl = [0; 0; 0; 0];                               
% ql = [0; 0; 0; 0];       
% pr = [0; 0; 0; 0];       
% %pr = [0; 0; 0; 0; ur(5)-dnt];
% %pr = [0; 0; 0; 0; ur(5)-d_max/2*(sin((t/delt_t)*pi/100)+2)];
% %pr = [0; 0; 0; 0; ur(5)-(3-rem(floor(t/(delt_t*30)),3))];%忘记了是哪种情况了
% qr = [0; 0; 0; 0];       
end