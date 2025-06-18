 function y=ODEsolver(y0)
    
h=1/24;
y=y0;

rk1=intraODEsys(y);
rk2=intraODEsys(y+0.5*h*rk1);
rk3=intraODEsys(y+0.5*h*rk2);
rk4=intraODEsys(y+h*rk3);

y=y+1/6*(rk1+2*rk2+2*rk3+rk4)*h;    
return