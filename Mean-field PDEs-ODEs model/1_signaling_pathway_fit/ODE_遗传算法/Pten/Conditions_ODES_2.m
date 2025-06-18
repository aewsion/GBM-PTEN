global EGF IGF1 EGFR_I IGF1R_I PTEN AKT_I GSK_I
ExpData;

EGF=1;
IGF1=1;
EGFR_I=0;
IGF1R_I=0;
PTEN=0;
AKT_I=0;
GSK_I=0;
x04=[pY1068_EGFR(1,1);p_ERK(1,1);p_IGF1Rbeta(1,1);1.5;0.5456;0.6624;0.5086];
Y1=zeros(24,7);%GSK IRF1 Gal9
AKT_I=0;
GSK_I=0;
t2=0:1/24:1;
[t2,y] = ode15s(@intraODEsys2,t2,x04);
Y1=y(25,:);%0,1/3,2/3,1,2,4,6,12,24



