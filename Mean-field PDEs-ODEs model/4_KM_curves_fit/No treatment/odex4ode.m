function dudt_ode = odex4ode(x,t,u)
global Ptest
global select_EGFR_I select_IGF1R_I dmax2 dmax3 delt_t TimeLength PTEN AKT_I GSK_I
V1=Ptest(1);
V2=Ptest(2);
V21=Ptest(3);
V22=Ptest(4);
V3=Ptest(5);
V4=Ptest(6);
V5=Ptest(7);
V51=Ptest(8);
V6=Ptest(9);
V7=Ptest(10);
d1=Ptest(11);
d2=Ptest(12);
d3=Ptest(13);
d4=Ptest(14);
d5=Ptest(15);
d6=Ptest(16);
d7=Ptest(17);
K11=Ptest(18);
K43=Ptest(19);
K21=Ptest(20);
K22=Ptest(21);
K51=Ptest(22);
K31=Ptest(23);
K52=Ptest(24);
K71=Ptest(25);
K42=Ptest(26);
K12=Ptest(27);
K23=Ptest(28);
K32=Ptest(29);
K41=Ptest(30);
K61=Ptest(31);
K13=Ptest(32);
K33=Ptest(33);
K44=Ptest(34);
n1=floor(Ptest(35));
n2=floor(Ptest(36));
% n2=1;
n3=floor(Ptest(37));
% n3=1;
n4=floor(Ptest(38));
n5=floor(Ptest(39));
n6=floor(Ptest(40));
n7=floor(Ptest(41));
n8=floor(Ptest(42));
n9=floor(Ptest(43));
K53=Ptest(44);
K54=Ptest(45);
K55=Ptest(46);


nx=length(x);
EGFR_I=zeros(1,nx);
IGF1R_I=zeros(1,nx);
% for i=1:1:nx
%     EGFR_I(i)=Drug2(x(i),t,select_EGFR_I,dmax2); %function aaa = Drug2(x,t,select,dmax,dmax_1)
%     IGF1R_I(i)=Drug2(x(i),t,select_IGF1R_I,dmax3);
% end

up=zeros(7,nx);
for jj=1:1:nx
    up(1,jj)=V1*u(6,jj)/(K11+u(6,jj))/(1+u(9,jj)/K12)/(1+EGFR_I(jj)/K13)*(1-u(8,jj))-d1*u(8,jj);                   %EGFR激活的
    up(2,jj)=V2*(1+V21*u(8,jj)^n1/(K21^n1+u(8,jj)^n1))*(1+V22*u(10,jj)/(K22+u(10,jj)))/(1+u(11,jj)/K23)*(1-u(9,jj))-d2*u(9,jj);    %ERK激活的
    up(3,jj)=V3*u(7,jj)/(K31+u(7,jj))/(1+u(9,jj)/K32)/(1+IGF1R_I(jj)/K33)*(1-u(10,jj))-d3*u(10,jj);                %IGFR激活的
    up(4,jj)=V4/(1+PTEN/K41)*u(8,jj)/(K42+u(8,jj))*u(10,jj)/(K43+u(10,jj))/(1+AKT_I^n3/K44^n3)*(1-u(11,jj))-d4*u(11,jj);%AKT激活的
    up(5,jj)=V5*u(11,jj)^n2/(K51^n2+u(11,jj)^n2)*(1+V51*GSK_I^n7/(K52^n7+GSK_I^n7))*(1-u(12,jj))-d5*u(12,jj);%GSK
    up(6,jj)=V6*(1-u(13,jj))-d6*u(13,jj)/(1+u(12,jj)^n8/K61^n8);%IRF1
    up(7,jj)=V7*u(13,jj)^n9/(K55^n9+u(13,jj)^n9)*(1-u(14,jj))-d7*u(14,jj); %Gal9

end
dudt_ode=up;
end