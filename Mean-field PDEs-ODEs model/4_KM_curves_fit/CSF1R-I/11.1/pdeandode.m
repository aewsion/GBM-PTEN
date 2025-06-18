function dudt = pdeodes(tnow,y)
global V1 a V21 V22 V3 V4 K11 d1 d2 d3 d4 K12 K21 K22 K23 K31 K32 K41 K42 K13 K33 n
global EGFR_I IGF1R_I
global npde node nx singular m xmesh xi varargin xim zxmp1 xzmp1
  %PDEODES  Assemble the difference equations and evaluate the time derivative
  %   for the ODE system.

    u = reshape(y,npde+node,nx);%reshape(A,[2,3]) 将 A 重构为一个 2×3 矩阵
    u1=u(1:npde,1:nx);
    u2=u(npde+1:npde+node,1:nx);
    up = zeros(npde+node,nx);
    [U1,Ux1] = pdentrp(singular,m,xmesh(1),u(1:npde+node,1),xmesh(2),u(1:npde+node,2),xi(1));
    U=zeros(npde+node,1);
    Ux11=Ux1(1:npde,1);
    U=U1;
    [cL,fL,sL] = feval(pde,xi(1),tnow,U,Ux11,varargin{:});


    %  Evaluate the boundary conditions
    [pL,qL,pR,qR] = feval(bc,xmesh(1),u(:,1),xmesh(nx),u(:,nx),tnow,varargin{:});

    %  Left boundary
    if singular
      denom = cL;
      denom(denom == 0) = 1;
      up(1:npde,1) = (sL + (m+1) * fL / xi(1)) ./ denom;
    else
      up(1:npde,1) = pL;
      idx = (qL ~= 0);
      denom = (qL(idx)/xmesh(1)^m) .* (zxmp1(1)*cL(idx));
      denom(denom == 0) = 1;
      up(idx,1) = ( pL(idx) + (qL(idx)/xmesh(1)^m) .* ...
                    (xim(1)*fL(idx) + zxmp1(1)*sL(idx))) ./ denom;
    end
    %  Interior points
    for ii = 2:nx-1
      [U1,Ux1] = pdentrp(singular,m,xmesh(ii),u(1:npde+node,ii),xmesh(ii+1),u(1:npde+node,ii+1),xi(ii));
      U = U1;
      Ux11=Ux1(1:npde);
      [cR,fR,sR] = feval(pde,xi(ii),tnow,U,Ux11,varargin{:});
%       [U1,Ux1] = pdentrp(singular,m,xmesh(ii),u1(1:npde,ii),xmesh(ii+1),u1(1:npde,ii+1),xi(ii));
%       U(1:npde) = U1;
%       U(npde+1:npde+node) = u2(:,ii);
%       [cR,fR,sR] = feval(pde,xi(ii),tnow,U,Ux1,varargin{:});        
      denom = zxmp1(ii) * cR + xzmp1(ii-1) * cL;
      denom(denom == 0) = 1;
      up(1:npde,ii) = ((xim(ii) * fR - xim(ii-1) * fL) + ...
                  (zxmp1(ii) * sR + xzmp1(ii-1) * sL)) ./ denom;

      cL = cR;
      fL = fR;
      sL = sR;
    end
    %  Right boundary
    up(1:npde,nx) = pR;
    idx = (qR ~= 0);
    denom = -(qR(idx)/xmesh(nx)^m) .* (xzmp1(nx-1)*cL(idx));
    denom(denom == 0) = 1;
    up(idx,nx) = ( pR(idx) + (qR(idx)/xmesh(nx)^m) .* ...
                   (xim(nx-1)*fL(idx) - xzmp1(nx-1)*sL(idx))) ./ denom;
    for jj = 1:nx         
    up(npde+1,jj) = V1*u(5,jj)/5/(K11+u(5,jj)/5)/(1+u(npde+2,jj)/K12)/(1+EGFR_I/K13)*(1-u(npde+1,jj))-d1*u(npde+1,jj);
    up(npde+2,jj) = a*(1+V21*u(npde+1,jj)^n/(K21^n+u(npde+1,jj)^n))*(1+V22*u(npde+3,jj)/(K22+u(npde+3,jj)))/(1+u(npde+4,jj)/K23)*(1-u(npde+2,jj))-d2*u(npde+2,jj); 
    up(npde+3,jj) = V3*u(6,jj)/3000/(K31+u(6,jj)/3000)/(1+u(npde+2,jj)/K32)/(1+IGF1R_I/K33)*(1-u(npde+3,jj))-d3*u(npde+3,jj);
    up(npde+4,jj) = V4*u(npde+1,jj)/(K41+u(npde+1,jj))*u(npde+3,jj)/(K42+u(npde+3,jj))*(1-u(npde+4,jj))-d4*u(npde+4,jj); 
    end
               

    dudt = up(:);
  end  % pdeodes
    
% --------------------------------------------------------------------------

end  % pdepe



function dudt = pdeodes(tnow,y)
  %PDEODES  Assemble the difference equations and evaluate the time derivative
  %   for the ODE system.

    u = reshape(y,npde+npde2,nx);%reshape(A,[2,3]) 将 A 重构为一个 2×3 矩阵
%     u1=u(1:npde,1:nx);
%     u2=u(npde+1:npde+npde2,1:nx);
%     up = zeros(npde+npde2,nx);
%     [U1,Ux1] = pdentrp(singular,m,xmesh(1),u1(1:npde,1),xmesh(2),u1(1:npde,2),xi(1));
%     U=zeros(npde+npde2,1);
%     U(1:npde) = U1;
%     U(npde+1:npde+npde2,1) = u2(:,1);
%     [cL,fL,sL] = feval(pde,xi(1),tnow,U,Ux1,varargin{:});
    u1=u(1:npde,1:nx);
    u2=u(npde+1:npde+npde2,1:nx);
    up = zeros(npde+npde2,nx);
    [U1,Ux1] = pdentrp(singular,m,xmesh(1),u(1:npde+npde2,1),xmesh(2),u(1:npde+npde2,2),xi(1));
    U=zeros(npde+npde2,1);
    Ux11=Ux1(1:npde,1);
    U=U1;
    [cL,fL,sL] = feval(pde,xi(1),tnow,U,Ux11,varargin{:});


    %  Evaluate the boundary conditions
    [pL,qL,pR,qR] = feval(bc,xmesh(1),u(:,1),xmesh(nx),u(:,nx),tnow,varargin{:});

    %  Left boundary
    if singular
      denom = cL;
      denom(denom == 0) = 1;
      up(1:npde,1) = (sL + (m+1) * fL / xi(1)) ./ denom;
    else
      up(1:npde,1) = pL;
      idx = (qL ~= 0);
      denom = (qL(idx)/xmesh(1)^m) .* (zxmp1(1)*cL(idx));
      denom(denom == 0) = 1;
      up(idx,1) = ( pL(idx) + (qL(idx)/xmesh(1)^m) .* ...
                    (xim(1)*fL(idx) + zxmp1(1)*sL(idx))) ./ denom;
    end
    %  Interior points
    for ii = 2:nx-1
      [U1,Ux1] = pdentrp(singular,m,xmesh(ii),u(1:npde+npde2,ii),xmesh(ii+1),u(1:npde+npde2,ii+1),xi(ii));
      U = U1;
      Ux11=Ux1(1:npde);
      [cR,fR,sR] = feval(pde,xi(ii),tnow,U,Ux11,varargin{:});
%       [U1,Ux1] = pdentrp(singular,m,xmesh(ii),u1(1:npde,ii),xmesh(ii+1),u1(1:npde,ii+1),xi(ii));
%       U(1:npde) = U1;
%       U(npde+1:npde+npde2) = u2(:,ii);
%       [cR,fR,sR] = feval(pde,xi(ii),tnow,U,Ux1,varargin{:});        
      denom = zxmp1(ii) * cR + xzmp1(ii-1) * cL;
      denom(denom == 0) = 1;
      up(1:npde,ii) = ((xim(ii) * fR - xim(ii-1) * fL) + ...
                  (zxmp1(ii) * sR + xzmp1(ii-1) * sL)) ./ denom;

      cL = cR;
      fL = fR;
      sL = sR;
    end
    %  Right boundary
    up(1:npde,nx) = pR;
    idx = (qR ~= 0);
    denom = -(qR(idx)/xmesh(nx)^m) .* (xzmp1(nx-1)*cL(idx));
    denom(denom == 0) = 1;
    up(idx,nx) = ( pR(idx) + (qR(idx)/xmesh(nx)^m) .* ...
                   (xim(nx-1)*fL(idx) - xzmp1(nx-1)*sL(idx))) ./ denom;
    for jj = 1:nx         
    up(npde+1,jj) = V1*u(5,jj)/5/(K11+u(5,jj)/5)/(1+u(npde+2,jj)/K12)/(1+EGFR_I/K13)*(1-u(npde+1,jj))-d1*u(npde+1,jj);
    up(npde+2,jj) = a*(1+V21*u(npde+1,jj)^n/(K21^n+u(npde+1,jj)^n))*(1+V22*u(npde+3,jj)/(K22+u(npde+3,jj)))/(1+u(npde+4,jj)/K23)*(1-u(npde+2,jj))-d2*u(npde+2,jj); 
    up(npde+3,jj) = V3*u(6,jj)/3000/(K31+u(6,jj)/3000)/(1+u(npde+2,jj)/K32)/(1+IGF1R_I/K33)*(1-u(npde+3,jj))-d3*u(npde+3,jj);
    up(npde+4,jj) = V4*u(npde+1,jj)/(K41+u(npde+1,jj))*u(npde+3,jj)/(K42+u(npde+3,jj))*(1-u(npde+4,jj))-d4*u(npde+4,jj); 
    end
               

    dudt = up(:);
  end  % pdeodes

