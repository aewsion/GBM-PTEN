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
S1=Y1([1:5],[4,5,7]);
S2=Y2([1,2,3,5,9],[4,5,7]);
S3=Y1([6:14],[5,6,7]);
S4=Y3([1,2,3,4,7,13,19,37,73],[5,6,7]);


ObjV(iii)=sum(sum((Y4(:,1)-pY1068_EGFR).^2))+sum(sum((Y4(:,2)-p_ERK).^2))+sum(sum((Y4(:,3)-p_IGF1Rbeta).^2))+sum(sum((Y4(:,4)-pS473_AKT).^2))+sum(sum((S1-Exp1).^2))+sum(sum((S2-Exp2).^2))+sum(sum((S3-Exp3).^2))+sum(sum((S4-Exp4).^2))+(Y4(7,6)-0.2)^2+(Y4(7,5)-0.2)^2;   % y is simulated data; y0 is experimental data
end
% plot(ObjV)