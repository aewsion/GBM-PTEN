function u0 = pdex4ic(x)
global L
global C_T_0 C_M1_0 C_M2_0 CSF1_0 EGF_0 IGF1_0
global x_0

% ExpData;
% EGFR_0=pY1068_EGFR(1,1);%initial value
% ERK_0=p_ERK(1,1);
% IGF1R_0=p_IGF1Rbeta(1,1);
% AKT_0=pS473_AKT(1,1);

u0 = [(1-x/L)*0.5+0.1;
    x/L*0.05+0.12;
    x/L.*0.05+0.12;
    (1-x/L)*0.1+0.4;
    (1-x/L)*0.2+0.3;
    (1-x/L)*0.2+0.3;
    0;%IGF1
    0.7;%EGFR
    0.72;%ERK
    0;%IGF1R
    0.01;%AKT
    0.4;
    0.3;
    (1-x/L)*0.1+0.1%Gal9
    ];   %用药的时候  CSF1浓度改为了0.4

% u0 = [(1-x/L)*(0.1+0*0.1.*rand(1,length(x)))+0.5;
%     x/L*(0.05+0*0.1.*randn(1,length(x)))+0.12;
%     x/L.*(0.05+0*0.1.*randn(1,length(x)))+0.88;
%     (1-x/L)*(0.1+0*0.1.*rand(1,length(x)))+0.4;
%     (1-x/L)*(0.2+0*0.1.*rand(1,length(x)))+0.3;
%     (1-x/L)*(0.2+0*0.1.*rand(1,length(x)))+0.3;
%     0;%IGF1
%     0.7;%EGFR
%     0.72;%ERK
%     0;%IGF1R
%     0.01;%AKT
%     0.4;
%     0.3;
%     0.1
%     ];   %用药的时候  CSF1浓度改为了0.4

end