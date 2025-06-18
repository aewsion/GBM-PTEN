function [c,f1,s1] = pdex4pde(x,t,u,DuDx)
global A_T_EGF A_T_IGF1 A_M1 A_M2 epsilon r_T a_ERK K1 K2 K31 K32 K4 a_AKT  d_T d_TM1 a_M2M1 a_CSF1 a_A B_1 B_2 b_EGF
global CSF1_max EGF_max IGF1_max
global ERK AKT 
global x_0
global Ptest1
global delt_t
global select_CSF1R_I dmax1
global D_d1 eta1
% global drug
% Parameters_nondimensional_4_gai_jin_fangcheng_5
% c = [1; 1; 1; 1; 1; 1];  
c=ones(7,1);

% EGFR_I=0;
% IGF1R_I=0;
% x1=x/x_0;

% bbb=abs(u)<=0.005;
%     u(bbb)=0.005;
t
% f1 = [1,0,0,0,-A_T_EGF*u(1),-A_T_IGF1*u(1);%A_T是新的A1了
%       0,1,0,-A_M1*u(2),0,0;%A是新置换的巨噬细胞的参数
%       0,0,1,-A_M2*u(3),0,0;
%       0,0,0,1/epsilon,0,0;
%       0,0,0,0,1/epsilon,0;
%       0,0,0,0,0,1/epsilon]*DuDx; 
f1=[D_d1
    ]*DuDx;

% f = [1,0,0,0,-A_T_EGF*u(1),-A_T_IGF1*u(1),0,0,0,0;%A_T是新的A1了
%     0,1,0,-A_M1*u(2),0,0,0,0,0,0;%A是新置换的巨噬细胞的参数
%     0,0,1,-A_M2*u(3),0,0,0,0,0,0;
%     0,0,0,1/epsilon,0,0,0,0,0,0;
%     0,0,0,0,1/epsilon,0,0,0,0,0;
%     0,0,0,0,0,1/epsilon,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0]*DuDx; 

% f = [0,0,0,0;
%     0,0,0,0;
%     0,0,0,0;
%     0,0,0,0]*DuDx; 
% % 
%      drug=0;%边界条件药物为0
%      IL4=0; 
% if t>=delt_t*24*25
%      t1=t-delt_t*24*25;
%      [drug,buyao] = picture(x,t1); %边界条件一药物的值    
%      IL4=IL4_1(x,t1);
% else
%     drug=0;
%     IL4=0; 
% end
%      [drug,buyao] = picture(x,t); %边界条件一药物的值    
%      IL4=IL4_1(x,t);
%      [buyao,drug] = picture(x,t); %边界条件二的时候药物的值
%      IL4=IL4_2(x,t);
%      drug=Drug_34(x,t);%边界条件三的时候药物的值 
%      IL4=IL4_3(x,t);

% %调试
% %      R1=r_T*(1+a_r*u(5)./(K1+u(5))+a_r*u(6)./(K4+u(6)));  % growth rate
%      R1=r_T*(1+a_ERK*value_ERK/(K5+value_ERK)+a_AKT*value_AKT/(K6+value_AKT));
%      D1=d_T*(1+d_TM1*u(2));  % death rate    
% %      R3=a_M2M1*u(3)-u(2).*(a_CSF1*u(4)./(K21+K22*drug+u(4))+a_A*IL4/(K3+IL4));
%      R3=a_M2M1*u(3)-u(2)*(a_CSF1*u(4)/(K21+K22*drug+u(4))+a_A*IL4/(K3+IL4));
%      bs=u(3)*IL4-u(6);
%     if bs<0
%         bs=0;
%     end
%      drug_CSF1R_I=Drug2(x,t,select_CSF1R_I,dmax1);
%      IL4=A_Drugsimulation2(x,t,select_CSF1R_I,dmax1);
% %调试
% %      R1=r_T*(1+a_r*u(5)./(K1+u(5))+a_r*u(6)./(K4+u(6)));  % growth rate
%      R1=r_T*(1+a_ERK*u(8)/(K1+u(8))+a_AKT*u(10)/(K2+u(10)));
%      D1=d_T*(1+d_TM1*u(2));  % death rate    
% %      R3=a_M2M1*u(3)-u(2).*(a_CSF1*u(4)./(K21+K22*drug_CSF1R_I+u(4))+a_A*IL4/(K3+IL4));
%      R3=a_M2M1*u(3)-u(2)*(a_CSF1*u(4)/(K31+K32*drug_CSF1R_I+u(4))+a_A*IL4/(K4+IL4));
% %      bs=u(3)*IL4-u(6);
% %     if bs<0
% %         bs=0;
% %     end

%     R7=V1*EGF/(K11+EGF)/(1+u(2)/K12)/(1+EGFR_I/K13)*(1-u(1))-d1*u(1);
%     R8=a*(1+V21*u(1)^n/(K21^n+u(1)^n))*(1+V22*u(3)/(K22+u(3)))/(1+u(4)/K23)*(1-u(2))-d2*u(2);
%     R9=V3*IGF1/(K31+IGF1)/(1+u(2)/K32)/(1+IGF1R_I/K33)*(1-u(3))-d3*u(3);
%     R10=V4*u(1)/(K41+u(1))*u(3)/(K42+u(3))*(1-u(4))-d4*u(4);
% 新的
%  s = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)-u(4));1/epsilon*B_2*(u(3).*u(4)./(K21+K22*drug_CSF1R_I+u(4))-u(5));1/epsilon*bs;
%      R7;R8;R9;R10]; 
%  s1 = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)*(1-u(4)/CSF1_max)-u(4));1/epsilon*B_2*(b_EGF*u(3).*u(4)./(K210+K220*drug_CSF1R_I+u(4))*(1-u(5)/EGF_max)-u(5));1/epsilon*(u(3)*IL4*(1-u(6)/IGF1_max)-u(6))];
% s1 = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)*(1-u(4)/CSF1_max)-u(4));1/epsilon*B_2*(b_EGF*u(3)*(1-u(5)/EGF_max)-u(5));1/epsilon*(u(3)*IL4*(1-u(6)/IGF1_max)-u(6))]; 
% %  s2 = [R7;R8;R9;R10]; 
s1=[-u(1)*eta1];
end