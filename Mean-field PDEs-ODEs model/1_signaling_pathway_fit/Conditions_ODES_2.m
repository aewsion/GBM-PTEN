global EGF IGF1 EGFR_I IGF1R_I PTEN AKT_I GSK_I Ptest1
Conditions;
ExpData;
Ptest1;

x1=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);pS473_AKT(1,1);0.1;0.1;0.1];%Gal9 考虑要不要改
Y4=zeros(11,6);
for ii=1:1:11
    EGF=Cond1(ii,1);
    IGF1=Cond1(ii,2);
    EGFR_I=Cond1(ii,3);
    IGF1R_I=Cond1(ii,4);
    PTEN=1;
    AKT_I=0;
    GSK_I=0;
    t2=0:1/24:1;
    [t2,y] = ode15s(@intraODEsys2,t2,x1);
    Y4(ii,1)=y(25,1);
    Y4(ii,2)=y(25,2);
    Y4(ii,3)=y(25,3);
    Y4(ii,4)=y(25,4);
    Y4(ii,5)=y(25,5);
    Y4(ii,6)=y(25,7);
end
x011=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);2.1673;0.4955;0.6;0.6];
x01=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);2.1673;0.4955;0.4;0.6];
x02=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);1.5;0.5456;0.6624;0.5086];
Y1=zeros(14,7);
EGF=1;
IGF1=1;
EGFR_I=0;
IGF1R_I=0;
PTEN=0;
for ii=1:30
    AKT_I=Cond2(1,ii);
    GSK_I=Cond2(2,ii);
    t2=0:1/24:1;
    if ii==1
        [t2,y] = ode15s(@intraODEsys2,t2,x01);
    elseif (1<ii)&&(ii<=21)
        [t2,y] = ode15s(@intraODEsys2,t2,x01);
    elseif ii==22
        [t2,y] = ode15s(@intraODEsys2,t2,x02);
    else
        [t2,y] = ode15s(@intraODEsys2,t2,x02);
    end
    Y1(ii,:)=y(25,:);
end

x03=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);2.1673;0.4955;0.6;0.8517];

AKT_I=2;
GSK_I=0;
t2=0:1/80:2;
[t2,y] = ode15s(@intraODEsys2,t2,x03);
Y2=y; %0,6,12,24,48


% x04=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);0.7;0.58;0.74;0.39];
x04=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);1.5;0.5456;0.6624;0.5086];
Y3=zeros(73,7);%GSK IRF1 Gal9
AKT_I=0;
GSK_I=5;
t2=0:1/72:1;
[t2,y] = ode15s(@intraODEsys2,t2,x04);
Y3=y;%0,1/3,2/3,1,2,4,6,12,24



