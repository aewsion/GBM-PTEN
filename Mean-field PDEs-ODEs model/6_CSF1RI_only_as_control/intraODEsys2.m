function xprime=intraODEsys2(t,y)


global Ptest1
global EGF IGF1 EGFR_I IGF1R_I PTEN AKT_I GSK_I
%model 1
V1=Ptest1(1);
V2=Ptest1(2);
V21=Ptest1(3);
V22=Ptest1(4);
V3=Ptest1(5);
V4=Ptest1(6);
V5=Ptest1(7);
V51=Ptest1(8);
V6=Ptest1(9);
V7=Ptest1(10);
d1=Ptest1(11);
d2=Ptest1(12);
d3=Ptest1(13);
d4=Ptest1(14);
d5=Ptest1(15);
d6=Ptest1(16);
d7=Ptest1(17);
K11=Ptest1(18);
K43=Ptest1(19);
K21=Ptest1(20);
K22=Ptest1(21);
K51=Ptest1(22);
K31=Ptest1(23);
K52=Ptest1(24);
K71=Ptest1(25);
K42=Ptest1(26);
K12=Ptest1(27);
K23=Ptest1(28);
K32=Ptest1(29);
K41=Ptest1(30);
K61=Ptest1(31);
K13=Ptest1(32);
K33=Ptest1(33);
K44=Ptest1(34);
n1=floor(Ptest1(35));
n2=floor(Ptest1(36));
% n2=1;
n3=floor(Ptest1(37));
% n3=1;
n4=floor(Ptest1(38));
n5=floor(Ptest1(39));
 n6=floor(Ptest1(40));
n7=floor(Ptest1(41));
n8=floor(Ptest1(42));
n9=floor(Ptest1(43));
K53=Ptest1(44);
K54=Ptest1(45);
K55=Ptest1(46);


% 新方程 model 1
xprime=[
    V1*EGF/(K11+EGF)/(1+y(2)/K12)/(1+EGFR_I/K13)*(1-y(1))-d1*y(1);                    %EGFR激活的
    V2*(1+V21*y(1)^n1/(K21^n1+y(1)^n1))*(1+V22*y(3)/(K22+y(3)))/(1+y(4)/K23)*(1-y(2))-d2*y(2);    %ERK激活的
    V3*IGF1/(K31+IGF1)/(1+y(2)/K32)/(1+IGF1R_I/K33)*(1-y(3))-d3*y(3);                %IGFR激活的
    V4/(1+PTEN/K41)*y(1)/(K42+y(1))*y(3)/(K43+y(3))/(1+AKT_I^n3/K44^n3)*(1-y(4))-d4*y(4);%AKT激活的
    V5*y(4)^n2/(K51^n2+y(4)^n2)*(1+V51*GSK_I^n7/(K52^n7+GSK_I^n7))*(1-y(5))-d5*y(5);%GSK
    V6*(1-y(6))-d6*y(6)/(1+y(5)^n8/K61^n8);%IRF1
    V7*y(6)^n9/(K55^n9+y(6)^n9)*(1-y(7))-d7*y(7); %Gal9
    ];



return