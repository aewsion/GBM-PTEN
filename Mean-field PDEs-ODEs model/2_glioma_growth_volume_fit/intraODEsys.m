function xprime=intraODEsys(y)


global Ptest1 
global EGF IGF1 EGFR_I IGF1R_I 
%model 1
V1=Ptest1(1);
a=Ptest1(2);
V21=Ptest1(3);
V22=Ptest1(4);
V3=Ptest1(5);
V4=Ptest1(6);
d1=Ptest1(7);
d2=Ptest1(8);
d3=Ptest1(9);
d4=Ptest1(10);
K11=Ptest1(11);
K12=Ptest1(12);
K21=Ptest1(13);
K22=Ptest1(14);
K23=Ptest1(15);
K31=Ptest1(16);
K32=Ptest1(17);
K41=Ptest1(18);
K42=Ptest1(19);
K13=Ptest1(20);
K33=Ptest1(21);
n=floor(Ptest1(22));
% n2=floor(Ptest1(25));
% n3=floor(Ptest1(26));

% %model 2
% V1=Ptest1(1);
% a=Ptest1(2);
% V21=Ptest1(3);
% V22=Ptest1(4);
% V3=Ptest1(5);
% V4=Ptest1(6);
% d1=Ptest1(7);
% d2=Ptest1(8);
% d3=Ptest1(9);
% d4=Ptest1(10);
% K11=Ptest1(11);
% K12=Ptest1(12);
% K21=Ptest1(13);
% K22=Ptest1(14);
% K31=Ptest1(15);
% K41=Ptest1(16);
% K42=Ptest1(17);
% K13=Ptest1(18);
% K32=Ptest1(19);
% n=floor(Ptest1(20));

% %model 3
% V1=Ptest1(1);
% a=Ptest1(2);
% V21=Ptest1(3);
% V22=Ptest1(4);
% V3=Ptest1(5);
% V4=Ptest1(6);
% d1=Ptest1(7);
% d2=Ptest1(8);
% d3=Ptest1(9);
% d4=Ptest1(10);
% K11=Ptest1(11);
% K21=Ptest1(12);
% K22=Ptest1(13);
% K31=Ptest1(14);
% K32=Ptest1(15);
% K41=Ptest1(16);
% K42=Ptest1(17);
% K12=Ptest1(18);
% K33=Ptest1(19);
% n=floor(Ptest1(20));

% %model 4
% V1=Ptest1(1);
% a=Ptest1(2);
% V21=Ptest1(3);
% V22=Ptest1(4);
% V3=Ptest1(5);
% V4=Ptest1(6);
% d1=Ptest1(7);
% d2=Ptest1(8);
% d3=Ptest1(9);
% d4=Ptest1(10);
% K11=Ptest1(11);
% K21=Ptest1(12);
% K22=Ptest1(13);
% K31=Ptest1(14);
% K41=Ptest1(15);
% K42=Ptest1(16);
% K12=Ptest1(17);
% K32=Ptest1(18);
% n=floor(Ptest1(19));


% 新方程 model 1
xprime=[
    V1*EGF/(K11+EGF)/(1+y(2)/K12)/(1+EGFR_I/K13)*(1-y(1))-d1*y(1);                    %EGFR激活的
    a*(1+V21*y(1)^n/(K21^n+y(1)^n))*(1+V22*y(3)/(K22+y(3)))/(1+y(4)/K23)*(1-y(2))-d2*y(2);    %ERK激活的
    V3*IGF1/(K31+IGF1)/(1+y(2)/K32)/(1+IGF1R_I/K33)*(1-y(3))-d3*y(3);                %IGFR激活的
    V4*y(1)/(K41+y(1))*y(3)/(K42+y(3))*(1-y(4))-d4*y(4);                                 %AKT激活的
    ];
% % 新方程 model 2
% xprime=[
%     V1*EGF/(K11+EGF)/(1+y(2)/K12)/(1+EGFR_I/K13)*(1-y(1))-d1*y(1);                    %EGFR激活的
%     a*(1+V21*y(1)^n/(K21^n+y(1)^n))*(1+V22*y(3)/(K22+y(3)))*(1-y(2))-d2*y(2);    %ERK激活的
%     V3*IGF1/(K31+IGF1)/(1+IGF1R_I/K32)*(1-y(3))-d3*y(3);                %IGFR激活的
%     V4*y(1)/(K41+y(1))*y(3)/(K42+y(3))*(1-y(4))-d4*y(4);                                 %AKT激活的
%     ];
% % 新方程 model 3
% xprime=[
%     V1*EGF/(K11+EGF)/(1+EGFR_I/K12)*(1-y(1))-d1*y(1);                    %EGFR激活的
%     a*(1+V21*y(1)^n/(K21^n+y(1)^n))*(1+V22*y(3)/(K22+y(3)))*(1-y(2))-d2*y(2);    %ERK激活的
%     V3*IGF1/(K31+IGF1)/(1+y(2)/K32)/(1+IGF1R_I/K33)*(1-y(3))-d3*y(3);                %IGFR激活的
%     V4*y(1)/(K41+y(1))*y(3)/(K42+y(3))*(1-y(4))-d4*y(4);                                 %AKT激活的
%     ];
% % 新方程 model 4
% xprime=[
%     V1*EGF/(K11+EGF)/(1+EGFR_I/K12)*(1-y(1))-d1*y(1);                    %EGFR激活的
%     a*(1+V21*y(1)^n/(K21^n+y(1)^n))*(1+V22*y(3)/(K22+y(3)))*(1-y(2))-d2*y(2);    %ERK激活的
%     V3*IGF1/(K31+IGF1)/(1+IGF1R_I/K32)*(1-y(3))-d3*y(3);                %IGFR激活的
%     V4*y(1)/(K41+y(1))*y(3)/(K42+y(3))*(1-y(4))-d4*y(4);                                 %AKT激活的
%     ];



return 