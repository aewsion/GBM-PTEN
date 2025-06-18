function ObjV=Objfun1(rvChrom)
[m,n]=size(rvChrom);
global time
global Ptest1     % Get the original true data 
ObjV=zeros(m,1);
for iii=1:m

Ptest1=rvChrom(iii,:);%每次只算一条染色体rvChrom(i,:)
% Conditions_ODES;
Conditions_ODES_2;
%%%%%
time=time+1;
% %新的
% ObjV(i)=sum(sum((Y(:,1)-pY1068_EGFR).^2))+sum(sum((Y(:,2)-p_ERK).^2))+sum(sum((Y(:,3)-p_IGF1Rbeta).^2))+sum(sum((Y(:,4)-pS473_AKT).^2))-sum(sum((Y(6,2)-p_ERK(6,1)).^2));   % y is simulated data; y0 is experimental data
%新的
ObjV(iii)=sum(sum((Y2(:,1)-pY1068_EGFR).^2))+sum(sum((Y2(:,2)-p_ERK).^2))+sum(sum((Y2(:,3)-p_IGF1Rbeta).^2))+sum(sum((Y2(:,4)-pS473_AKT).^2));   % y is simulated data; y0 is experimental data

end
% plot(ObjV)