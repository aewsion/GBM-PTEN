    global EGF IGF1 EGFR_I IGF1R_I     
    t2=0:1/24:1;
    Conditions;
    ExpData;
    c1=size(p_IGF1Rbeta);
    c2=c1(1); 
    yy0=zeros(4,1);
    yy0(1,1)=pY1068_EGFR(1,1);%x0 initial value  4个ODE的初始值
    yy0(2,1)=p_ERK(1,1);
    yy0(3,1)=p_IGF1Rbeta(1,1);
    yy0(4,1)=pS473_AKT(1,1);
    Y2=zeros(c2,4);
    for ii=1:1:c2 
        EGF=Cond1(ii,1);
        IGF1=Cond1(ii,2);
        EGFR_I=Cond1(ii,3);          
        IGF1R_I=Cond1(ii,4);  
%         options=odeset('reltol',1e-3);%精度函数不能少
        [t2,y] = ode15s(@intraODEsys2,t2,yy0(:));
        Y2(ii,1)=y(25,1);
        Y2(ii,2)=y(25,2);
        Y2(ii,3)=y(25,3);
        Y2(ii,4)=y(25,4);
%         Q(ii,:,:)=y;
    end
    
    
    
    
% %     Q=zeros(a2,49,4);
%     x0=ones(4,1);
% %      Q(ii,1,1)=pY1068_EGFR(1,1);%initial value
% %      Q(ii,1,2)=p_ERK(1,1);
% %      Q(ii,1,3)=p_IGF1Rbeta(1,1);
% %      Q(ii,1,4)=pS473_AKT(1,1);
%      %下面是模拟两天的情况
%      for i=1:1:c2   
%          x0(1,1)=pY1068_EGFR(1,1);%x0 initial value
%          x0(2,1)=p_ERK(1,1);
%          x0(3,1)=p_IGF1Rbeta(1,1);
%          x0(4,1)=pS473_AKT(1,1);
%          EGF=Cond1(i,1);
%          IGF1=Cond1(i,2);
%          EGFR_I=Cond1(i,3);          
%          IGF1R_I=Cond1(i,4);  
%          y0=x0;
%          for t1=1:1:24
%          y=ODEsolver(y0);
%          y0=y;
%          end
%          Y(i,1)=y(1,1);
%          Y(i,2)=y(2,1);
%          Y(i,3)=y(3,1);
%          Y(i,4)=y(4,1);
%      end
% %  end    
