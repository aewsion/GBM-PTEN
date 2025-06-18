function u0 = pdex4ic(x)
global L
global C_T_0 C_M1_0 C_M2_0 CSF1_0 EGF_0 IGF1_0
global x_0
% x1=x*x_0;
% a=floor(x/x_0);
% b=x-x_0*a;
% a2=a+1;
% if b~=0
%     value1=C_T_0(1,a2)+b/x_0*(C_T_0(1,a2+1)-C_T_0(1,a2));
%     value2=C_M1_0(1,a2)+b/x_0*(C_M1_0(1,a2+1)-C_M1_0(1,a2));
%     value3=C_M2_0(1,a2)+b/x_0*(C_M2_0(1,a2+1)-C_M2_0(1,a2));
%     value4=CSF1_0(1,a2)+b/x_0*(CSF1_0(1,a2+1)-CSF1_0(1,a2));
%     value5=EGF_0(1,a2)+b/x_0*(EGF_0(1,a2+1)-EGF_0(1,a2));
%     value6=IGF1_0(1,a2)+b/x_0*(IGF1_0(1,a2+1)-IGF1_0(1,a2));
% else
%     value1=C_T_0(1,a2);
%     value2=C_M1_0(1,a2);
%     value3=C_M2_0(1,a2);
%     value4=CSF1_0(1,a2);
%     value5=EGF_0(1,a2);
%     value6=IGF1_0(1,a2);
% end
% ExpData;
%      EGFR_0=pY1068_EGFR(1,1);%initial value
%      ERK_0=p_ERK(1,1);
%      IGF1R_0=p_IGF1Rbeta(1,1);
%      AKT_0=pS473_AKT(1,1);
% u0 = [0.5.*Heaviside(11-x).*(1+0.1.*randn(1,length(x))); 100/100.*(1+0.1.*randn(1,length(x)))+0.0; 1.*(1+0.01.*randn(1,length(x))); 1.*(1+0.01.*randn(1,length(x))); d_max*x/100];    

% u0 = [(1-x1/100)*(0.1+0*0.1.*rand(1,length(x1)))+0.5;x1/100*(0.05+0*0.1.*randn(1,length(x1)))+0.12; x1/100.*(0.05+0*0.1.*randn(1,length(x1)))+0.88; (1-x1/100)*(0.1+0*0.1.*rand(1,length(x1)))+0.7; (1-x1/100)*(0.2+0*0.1.*rand(1,length(x1)))+0.3;0;
%     0.7;0.72;0.01;0.01];   %用药的时候
% u0 = [(1-x1/100)*(0.1+0*0.1.*rand(1,length(x1)))+0.5;x1/100*(0.05+0*0.1.*randn(1,length(x1)))+0.12; x1/100.*(0.05+0*0.1.*randn(1,length(x1)))+0.88; (1-x1/100)*(0.1+0*0.1.*rand(1,length(x1)))+0.4; (1-x1/100)*(0.2+0*0.1.*rand(1,length(x1)))+0.3;0;
%     0.7;0.72;0.01;0.01];   %用药的时候  CSF1浓度改为了0.4

u0=[x/L*1];
u0=[x1/100*1];

% u0 = [(1-x1/100)*(0.1)+0.08;x1/100*(0.05)+0.12; x1/100.*(0.05)+0.88; (1-x1/100)*(0.1+0*0.1.*rand(1,length(x1)))+0.7; (1-x1/100)*(0.2+0*0.1.*rand(1,length(x1)))+0.3;0;
%     0.7;0.72;0.01;0.01];   %从小开始生长


% u0 = [0.1; x/100.*(1+0.1.*randn(1,length(x)))+1; x/100.*(1+0.1.*randn(1,length(x)))+1;0.1; 0.1; 0;EGFR_0;ERK_0;IGF1R_0;AKT_0]; 
% u0 = [0.1; 1.5; 1.5;0.1; 0.1; 0;EGFR_0;ERK_0;IGF1R_0;AKT_0]; 
% u0 = [EGFR_0;ERK_0;IGF1R_0;AKT_0]; 
end