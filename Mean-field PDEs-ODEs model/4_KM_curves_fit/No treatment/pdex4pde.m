function [c,f1,s1] = pdex4pde(x,t,u,DuDx)
global A_EGF A_IGF1 A_CSF1 A_Gal9 epsilon r_T alpha_ERK alpha_AKT d_T d_TM1 B_CSF1 B_Gal9 B_M1M2 
global K1 K2 K3 K4 K5 b_EGF M1 M2 M3
global CSF1_max EGF_max IGF1_max
global select_CSF1R_I dmax1 select_Gal9_I dmax2  K_CSF1RI K_Gal9I K_A B_CSF1RI
% global drug
% Parameters_nondimensional_4_gai_jin_fangcheng_4
c = [1;1;1;1;1;1;1];%7pde 7ode
f1 = [1,0,0,0,0,-A_EGF*u(1),-A_IGF1*u(1);%A_T是新的A_EGF了
    0,1,0,-A_CSF1*u(2),-A_Gal9*u(2),0,0;%A是新置换的巨噬细胞的参数
    0,0,1,-A_CSF1*u(3),-A_Gal9*u(3),0,0;
    0,0,0,1/epsilon,0,0,0;
    0,0,0,0,1/epsilon,0,0;
    0,0,0,0,0,1/epsilon,0;
    0,0,0,0,0,0,1/epsilon]*DuDx;


drug_CSF1R_I=Drug2(x,t,select_CSF1R_I,dmax1);%返回药物
IL4=A_Drugsimulation2(x,t,select_CSF1R_I,dmax1);%药物的累积

drug_Gal9_I=Drug2(x,t,select_Gal9_I,dmax2);


IL4;
R1=r_T*(0+alpha_ERK/2*u(9)/(K4+u(9))+alpha_AKT*u(11)/(K5+u(11)));
D1=d_T*(1+d_TM1*u(2)); % death rate
% R3=B_M1M2*u(3)-u(2)*(B_CSF1*u(4)/(K1+u(4)+10*K_CSF1RI*drug_CSF1R_I)+B_Gal9*u(5)/(K2+u(5)+10*K_Gal9I*drug_Gal9_I)+B_CSF1RI*IL4/(K_A/10+IL4));%M1 M2之间的转换
% R3=B_M1M2*u(3)-u(2)*(B_CSF1*u(4)/(K1+u(4)+K_CSF1RI*drug_CSF1R_I)+B_Gal9*u(5)/(K2+u(5)+10*K_Gal9I*drug_Gal9_I)+IL4);%M1 M2之间的转换
u(3);
R3=B_M1M2*u(3)-u(2)*(B_CSF1*u(4)/(K1+u(4)+K_CSF1RI*drug_CSF1R_I/2)*10*1.5+10/1.5*B_Gal9*u(5)/(K2+u(5)+10*K_Gal9I*drug_Gal9_I)+500*IL4/(1+IL4));
% a1=u(2)
% a2=u(3)
% a3=u(5)
% a4=u(11)


% M2
s1 = [R1*u(1)*(1-u(1))-D1*u(1);
    R3;
    -R3;
    M1/epsilon*(u(1)*(1-u(4)/CSF1_max)-u(4));
    M2/epsilon*(u(1)*u(14)/(K3+u(14))-u(5));
    M3/epsilon*(b_EGF*u(3)*(1-u(6)/EGF_max)-u(6));
    1/epsilon*(0.1*u(3)*(1-u(7)/IGF1_max)-u(7))];%C_T M1 M2 CSF1 EGF IGF1


end