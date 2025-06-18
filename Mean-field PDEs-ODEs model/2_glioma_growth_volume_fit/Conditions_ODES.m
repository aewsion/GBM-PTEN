    global EGF IGF1 EGFR_I IGF1R_I 
    Conditions;
    ExpData;
    c1=size(p_IGF1Rbeta);
    c2=c1(1);
    Y=zeros(c2,4);
%     Q=zeros(a2,49,4);
    x0=ones(4,1);
%      Q(ii,1,1)=pY1068_EGFR(1,1);%initial value
%      Q(ii,1,2)=p_ERK(1,1);
%      Q(ii,1,3)=p_IGF1Rbeta(1,1);
%      Q(ii,1,4)=pS473_AKT(1,1);
     %下面是模拟两天的情况
     for i=1:1:c2   
         x0(1,1)=pY1068_EGFR(1,1);%x0 initial value
         x0(2,1)=p_ERK(1,1);
         x0(3,1)=p_IGF1Rbeta(1,1);
         x0(4,1)=pS473_AKT(1,1);
         EGF=Cond1(i,1);
         IGF1=Cond1(i,2);
         EGFR_I=Cond1(i,3);          
         IGF1R_I=Cond1(i,4);  
         y0=x0;
         for t1=1:1:24
         y=ODEsolver(y0);
         y0=y;
         end
         Y(i,1)=y(1,1);
         Y(i,2)=y(2,1);
         Y(i,3)=y(3,1);
         Y(i,4)=y(4,1);
     end
%  end    
