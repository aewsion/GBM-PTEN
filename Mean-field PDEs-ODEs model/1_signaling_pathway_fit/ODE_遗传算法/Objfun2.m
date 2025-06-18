function ObjV=Objfun2(rvChrom)
global dmax1 dmax2 dmax3 period inter
[m,n]=size(rvChrom);
% global time
% global Ptest1     % Get the original true data 
ObjV=zeros(m,1);
for iii=1:m
dmax1=rvChrom(iii,1);
dmax2=rvChrom(iii,2);
% dmax3=rvChrom(iii,3);
% inter=floor(rvChrom(iii,3));
% period=floor(rvChrom(iii,3))+floor(rvChrom(iii,4));
inter=1;
period=10;
% Conditions_ODES;
Result_model_with_drug;
if length(sol)<24*TimeLength+1
    C_T_mean_end=100;
else
    u1 = sol(:,:,1);
    C_T_mean=sum(u1(:,:),2)/101; 
    C_T_mean_end=C_T_mean(24*TimeLength+1,1);
end
%%%%%
% time=time+1;
% %新的
% ObjV(i)=sum(sum((Y(:,1)-pY1068_EGFR).^2))+sum(sum((Y(:,2)-p_ERK).^2))+sum(sum((Y(:,3)-p_IGF1Rbeta).^2))+sum(sum((Y(:,4)-pS473_AKT).^2))-sum(sum((Y(6,2)-p_ERK(6,1)).^2));   % y is simulated data; y0 is experimental data
%新的
if C_T_mean_end<0
    C_T_mean_end=100;
end
ObjV(iii)=C_T_mean_end*TimeLength/2+dmax1*TimeLength+dmax2*TimeLength/2;   % y is simulated data; y0 is experimental data
end
% plot(ObjV)