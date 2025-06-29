function varargout = ode15s(ode,tspan,y0,options,varargin)
%ODE15S Solve stiff differential equations and DAEs, variable order method.
%   [TOUT,YOUT] = ODE15S(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates 
%   the system of differential equations y' = f(t,y) from time T0 to TFINAL 
%   with initial conditions Y0. ODEFUN is a function handle. For a scalar T 
%   and a vector Y, ODEFUN(T,Y) must return a column vector corresponding 
%   to f(t,y). Each row in the solution array YOUT corresponds to a time
%   returned in the column vector TOUT.  To obtain solutions at specific
%   times T0,T1,...,TFINAL (all increasing or all decreasing), use TSPAN =
%   [T0 T1 ... TFINAL].     
%   
%   [TOUT,YOUT] = ODE15S(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the ODESET function. See ODESET for details. Commonly used options
%   are scalar relative error tolerance 'RelTol' (1e-3 by default) and vector
%   of absolute error tolerances 'AbsTol' (all components 1e-6 by default).  
%   If certain components of the solution must be non-negative, use
%   ODESET to set the 'NonNegative' property to the indices of these
%   components.  The 'NonNegative' property is ignored for problems
%   where there is a  mass matrix. 
%   
%   The Jacobian matrix df/dy is critical to reliability and efficiency. Use
%   ODESET to set 'Jacobian' to a function handle FJAC if FJAC(T,Y) returns 
%   the Jacobian df/dy or to the matrix df/dy if the Jacobian is constant. 
%   If the 'Jacobian' option is not set (the default), df/dy is approximated 
%   by finite differences. Set 'Vectorized' 'on' if the ODE function is coded 
%   so that ODEFUN(T,[Y1 Y2 ...]) returns [ODEFUN(T,Y1) ODEFUN(T,Y2) ...]. 
%   If df/dy is a sparse matrix, set 'JPattern' to the sparsity pattern of
%   df/dy, i.e., a sparse matrix S with S(i,j) = 1 if component i of f(t,y)
%   depends on component j of y, and 0 otherwise.    
%
%   ODE15S can solve problems M(t,y)*y' = f(t,y) with mass matrix M(t,y). Use
%   ODESET to set the 'Mass' property to a function handle MASS if MASS(T,Y) 
%   returns the value of the mass matrix. If the mass matrix is constant, 
%   the matrix can be used as the value of the 'Mass' option. Problems with
%   state-dependent mass matrices are more difficult. If the mass matrix does
%   not depend on the state variable Y and the function MASS is to be called
%   with one input argument T, set 'MStateDependence' to 'none'. If the mass
%   matrix depends weakly on Y, set 'MStateDependence' to 'weak' (the
%   default) and otherwise, to 'strong'. In either case the function MASS is
%   to be called with the two arguments (T,Y). If there are many differential
%   equations, it is important to exploit sparsity: Return a sparse
%   M(t,y). Either supply the sparsity pattern of df/dy using the 'JPattern'
%   property or a sparse df/dy using the Jacobian property. For strongly
%   state-dependent M(t,y), set 'MvPattern' to a sparse matrix S with S(i,j)
%   = 1 if for any k, the (i,k) component of M(t,y) depends on component j of
%   y, and 0 otherwise.    
%
%   If the mass matrix is non-singular, the solution of the problem is
%   straightforward. See examples FEM1ODE, FEM2ODE, BATONODE, or
%   BURGERSODE. If M(t0,y0) is singular, the problem is a differential-
%   algebraic equation (DAE). ODE15S solves DAEs of index 1. DAEs have
%   solutions only when y0 is consistent, i.e., there is a yp0 such that
%   M(t0,y0)*yp0 = f(t0,y0). Use ODESET to set 'MassSingular' to 'yes', 'no',
%   or 'maybe'. The default of 'maybe' causes ODE15S to test whether M(t0,y0)
%   is singular. You can provide yp0 as the value of the 'InitialSlope'
%   property. The default is the zero vector. If y0 and yp0 are not
%   consistent, ODE15S treats them as guesses, tries to compute consistent
%   values close to the guesses, and then goes on to solve the problem. See
%   examples HB1DAE or AMP1DAE.  
%
%   [TOUT,YOUT,TE,YE,IE] = ODE15S(ODEFUN,TSPAN,Y0,OPTIONS) with the 'Events'
%   property in OPTIONS set to a function handle EVENTS, solves as above 
%   while also finding where functions of (T,Y), called event functions, 
%   are zero. For each function you specify whether the integration is 
%   to terminate at a zero and whether the direction of the zero crossing 
%   matters. These are the three column vectors returned by EVENTS: 
%   [VALUE,ISTERMINAL,DIRECTION] = EVENTS(T,Y). For the I-th event function: 
%   VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the integration 
%   is to terminate at a zero of this event function and 0 otherwise. 
%   DIRECTION(I)=0 if all zeros are to be computed (the default), +1 if only 
%   zeros where the event function is increasing, and -1 if only zeros where 
%   the event function is decreasing. Output TE is a column vector of times 
%   at which events occur. Rows of YE are the corresponding solutions, and 
%   indices in vector IE specify which event occurred.    
%   
%   SOL = ODE15S(ODEFUN,[T0 TFINAL],Y0...) returns a structure that can be
%   used with DEVAL to evaluate the solution or its first derivative at 
%   any point between T0 and TFINAL. The steps chosen by ODE15S are returned 
%   in a row vector SOL.x.  For each I, the column SOL.y(:,I) contains 
%   the solution at SOL.x(I). If events were detected, SOL.xe is a row vector 
%   of points at which events occurred. Columns of SOL.ye are the corresponding 
%   solutions, and indices in vector SOL.ie specify which event occurred. 
%
%   Example
%         [t,y]=ode15s(@vdp1000,[0 3000],[2 0]);   
%         plot(t,y(:,1));
%     solves the system y' = vdp1000(t,y), using the default relative error
%     tolerance 1e-3 and the default absolute tolerance of 1e-6 for each
%     component, and plots the first component of the solution.
%
%   See also ODE23S, ODE23T, ODE23TB, ODE45, ODE23, ODE113, ODE15I,
%            ODESET, ODEPLOT, ODEPHAS2, ODEPHAS3, ODEPRINT, DEVAL,
%            ODEEXAMPLES, VDPODE, BRUSSODE, HB1DAE, FUNCTION_HANDLE.

%   ODE15S is a quasi-constant step size implementation in terms of backward
%   differences of the Klopfenstein-Shampine family of Numerical
%   Differentiation Formulas of orders 1-5. The natural "free" interpolants
%   are used. Local extrapolation is not done. By default, Jacobians are
%   generated numerically.  

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997, and in
%   Solving Index-1 DAEs in MATLAB and Simulink, L. F. Shampine,
%   M. W. Reichelt, and J. A. Kierzenka, SIAM Review, 41-3, 1999. 

%   Mark W. Reichelt, Lawrence F. Shampine, and Jacek Kierzenka, 12-18-97
%   Copyright 1984-2021 The MathWorks, Inc.

solver_name = 'ode15s';

import matlab.internal.math.nowarn.mldivide

if nargin < 4
  options = [];
  if nargin < 3
    y0 = [];
    if nargin < 2
      tspan = [];
      if nargin < 1
        error(message('MATLAB:ode15s:NotEnoughInputs'));
      end  
    end
  end
end

% Stats
nsteps   = 0;
nfailed  = 0;
nfevals  = 0; 
npds     = 0;
ndecomps = 0;
nsolves  = 0;

[ode, odeIsFuncHandle, odeTreatAsMFile] = packageAsFuncHandle(ode);

% Output
output_sol = (odeIsFuncHandle && (nargout==1));      % sol = odeXX(...)
output_ty  = (~output_sol && (nargout > 0));  % [t,y,...] = odeXX(...)
% There might be no output requested...

sol = []; kvec = []; dif3d = []; 
if output_sol
  sol.solver = solver_name;
  sol.extdata.odefun = ode;
  sol.extdata.options = options;                       
  sol.extdata.varargin = varargin;  
end  

% Handle solver arguments

[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, ...
 options, threshold, rtol, normcontrol, normy, hmax, htry, htspan] = ...
    odearguments(odeIsFuncHandle, odeTreatAsMFile, solver_name, ode, tspan, y0, options, varargin);
nfevals = nfevals + 1;
one2neq = (1:neq);

% Handle the output
if nargout > 0
  outputFcn = odeget(options,'OutputFcn',[],'fast');
else
  outputFcn = odeget(options,'OutputFcn',@odeplot,'fast');
end
outputArgs = {};      
if isempty(outputFcn)
  haveOutputFcn = false;
else
  haveOutputFcn = true;
  outputs = odeget(options,'OutputSel',1:neq,'fast');
  if isa(outputFcn,'function_handle')  
    % With MATLAB 6 syntax pass additional input arguments to outputFcn.
    outputArgs = varargin;
  end  
end
refine = max(1,odeget(options,'Refine',1,'fast'));
if ntspan > 2
  outputAt = 'RequestedPoints';         % output only at tspan points
elseif refine <= 1
  outputAt = 'SolverSteps';             % computed points, no refinement
else
  outputAt = 'RefinedSteps';            % computed points, with refinement
  S = (1:refine-1) / refine;
end
printstats = strcmp(odeget(options,'Stats','off','fast'),'on');

% Handle the event function 
[haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
    odeevents(odeIsFuncHandle,ode,t0,y0,options,varargin);

% Handle the mass matrix
[Mtype, Mt, Mfun, Margs, dMoptions] = odemass(odeIsFuncHandle,ode,t0,y0,...
                                              options,varargin);

% Non-negative solution components
idxNonNegative = odeget(options,'NonNegative',[],'fast');
nonNegative = ~isempty(idxNonNegative);
if nonNegative  
  if Mtype == 0
    % Explicit ODE -- modify the derivative function
    [ode,thresholdNonNegative] = odenonnegative(ode,y0,threshold,idxNonNegative);
    f0 = ode(t0,y0,odeArgs{:});
    nfevals = nfevals + 1;
  else
    % Linearly implicit ODE/DAE -- ignore non-negativity constraints
    warning(message('MATLAB:ode15s:NonNegativeIgnoredForLinearlyImplicitSystems'));   
    nonNegative = false;
    idxNonNegative = [];
  end  
end

% Handle the Jacobian
[Jconstant,Jac,Jargs,Joptions] = ...
    odejacobian(odeIsFuncHandle,ode,t0,y0,options,varargin);
Janalytic = isempty(Joptions);
    
t = t0;
y = y0;

yp0_OK = false;
DAE = false;
RowScale = [];
if Mtype > 0
  nz = nnz(Mt);
  if nz == 0
    error(message('MATLAB:ode15s:MassMatrixAllZero'))
  end
   
  Msingular = odeget(options,'MassSingular','maybe','fast');
  switch Msingular
    case 'no',     DAE = false;
    case 'yes',    DAE = true;
    case 'maybe',  DAE = (eps*nz*condest(Mt) > 1);       
  end
   
  if DAE
    yp0 = odeget(options,'InitialSlope',[],'fast');
    if isempty(yp0)
      yp0_OK = false;
      yp0 = zeros(neq,1);  
    else
      yp0 = yp0(:);
      if length(yp0) ~= neq
        error(message('MATLAB:ode15s:YoYPoLengthMismatch'));
      end
      % Test if (y0,yp0) are consistent enough to accept.
      yp0_OK = (norm(Mt*yp0 - f0) <= 1e-3*rtol*max(norm(Mt*yp0),norm(f0)));
    end   
    if ~yp0_OK           % Must compute ICs, so classify them.
      if Mtype >= 3  % state dependent
        ICtype = 3;
      else  % M, M(t)
        % Test for a diagonal mass matrix.
        [r,c] = find(Mt);
        if isequal(r,c)   % diagonal
          ICtype = 1;
        elseif ~issparse(Mt) % not diagonal but full
          ICtype = 2;
        else  % sparse, not diagonal
          ICtype = 3;
        end
      end      
    end
  end
end
Mcurrent = true;
Mtnew = Mt;

% if not set via 'options', initialize constant Jacobian here
if Jconstant 
  if isempty(Jac) % use odenumjac
    [Jac,Joptions.fac,nF] = odenumjac(ode, {t0,y0,odeArgs{:}}, f0, Joptions); %#ok<CCAT>
    nfevals = nfevals + nF;
    npds = npds + 1;
  elseif ~isa(Jac,'numeric')  % not been set via 'options'  
    Jac = Jac(t0,y0,Jargs{:}); % replace by its value
    npds = npds + 1;
  end
end

maxk = odeget(options,'MaxOrder',5,'fast');
bdf = strcmp(odeget(options,'BDF','off','fast'),'on');

% Initialize method parameters.
G = [1; 3/2; 11/6; 25/12; 137/60];
if bdf
  alpha = [0; 0; 0; 0; 0];
else
  alpha = [-37/200; -1/9; -0.0823; -0.0415; 0];
end
invGa = 1 ./ (G .* (1 - alpha));
erconst = alpha .* G + (1 ./ (2:6)');
difU = [ -1, -2, -3, -4,  -5;           % difU is its own inverse!
          0,  1,  3,  6,  10;
          0,  0, -1, -4, -10;
          0,  0,  0,  1,   5;
          0,  0,  0,  0,  -1 ];
maxK = 1:maxk;
[kJ,kI] = meshgrid(maxK,maxK);
difU = difU(maxK,maxK);
maxit = 4;

% Get the initial slope yp. For DAEs the default is to compute
% consistent initial conditions.
if DAE && ~yp0_OK
  if ICtype < 3
    [y,yp,f0,dfdy,nFE,nPD,Jfac] = daeic12(ode,odeArgs,t,ICtype,Mt,y,yp0,f0,...
                                          rtol,Jconstant,Jac,Jargs,Joptions); 
  else    
    [y,yp,f0,dfdy,nFE,nPD,Jfac,dMfac] = daeic3(ode,odeArgs,tspan,htry,Mtype,Mt,Mfun,...
                                               Margs,dMoptions,y,yp0,f0,rtol,Jconstant,...
                                               Jac,Jargs,Joptions);   
    if ~isempty(dMoptions)
      dMoptions.fac = dMfac;
    end        
  end  
  if ~isempty(Joptions)
    Joptions.fac = Jfac;
  end    
  nfevals = nfevals + nFE;
  npds = npds + nPD;
  if Mtype >= 3
    Mt = Mfun(t,y,Margs{:});
    Mtnew = Mt;
    Mcurrent = true;
  end
else
  if Mtype == 0 
    yp = f0;
  elseif DAE && yp0_OK
    yp = yp0;
  else  
    % use overloaded subfunction mldivide which throws no warning.
    if issparse(Mt)
      [L,U,P,Q,R] = lu(Mt);            
      yp = Q * (U \ (L \ (P * (R \ f0))));      
    else
      [L,U,p] = lu(Mt,'vector');      
      yp = U \ (L \ f0(p));
    end  
    ndecomps = ndecomps + 1;              
    nsolves = nsolves + 1;                
  end
    
  if Jconstant
    dfdy = Jac;
  elseif Janalytic
    dfdy = Jac(t,y,Jargs{:});     
    npds = npds + 1;                            
  else   % Joptions not empty
    [dfdy,Joptions.fac,nF] = odenumjac(ode, {t,y,odeArgs{:}}, f0, Joptions); %#ok<CCAT>
    nfevals = nfevals + nF;    
    npds = npds + 1;                            
  end     
end
Jcurrent = true;

% hmin is a small number such that t + hmin is clearly different from t in
% the working precision, but with this definition, it is 0 if t = 0.
hmin = 16*eps*abs(t);

if isempty(htry)
  % Compute an initial step size h using yp = y'(t).
  if normcontrol
    wt = max(normy,threshold);
    rh = 1.25 * (norm(yp) / wt) / sqrt(rtol);  % 1.25 = 1 / 0.8
  else
    wt = max(abs(y),threshold);
    rh = 1.25 * norm(yp ./ wt,inf) / sqrt(rtol);
  end
  absh = min(hmax, htspan);
  if absh * rh > 1
    absh = 1 / rh;
  end
  absh = max(absh, hmin);
  
  if ~DAE
    % The error of BDF1 is 0.5*h^2*y''(t), so we can determine the optimal h.
    h = tdir * absh;
    tdel = (t + tdir*min(sqrt(eps)*max(abs(t),abs(t+h)),absh)) - t;
    f1 = ode(t+tdel,y,odeArgs{:});
    nfevals = nfevals + 1;                
    dfdt = (f1 - f0) ./ tdel;
    DfDt = dfdt + dfdy*yp;
    if normcontrol
      if Mtype > 0 
          if issparse(Mt)  
              rh = 1.25 * sqrt(0.5 * (norm(U \ (L \ (P * (R \ DfDt)))) / wt) / rtol);
          else
              rh = 1.25 * sqrt(0.5 * (norm(U \ (L \ DfDt(p))) / wt) / rtol);
          end
      else
        rh = 1.25 * sqrt(0.5 * (norm(DfDt) / wt) / rtol);
      end
    else
      if Mtype > 0
        if issparse(Mt)
          rh = 1.25*sqrt(0.5*norm((Q * (U \ (L \ (P * (R \ DfDt))))) ./ wt,inf) / rtol);
        else  
          rh = 1.25*sqrt(0.5*norm((U \ (L \ DfDt(p))) ./ wt,inf) / rtol);
        end  
      else
        rh = 1.25 * sqrt(0.5 * norm( DfDt ./ wt,inf) / rtol);
      end
    end
    absh = min(hmax, htspan);
    if absh * rh > 1
      absh = 1 / rh;
    end
    absh = max(absh, hmin);
  end
else
  absh = min(hmax, max(hmin, htry));
end
h = tdir * absh;

% Initialize.
k = 1;                                  % start at order 1 with BDF1
K = 1;                                  % K = 1:k
klast = k;
abshlast = absh;

dif = zeros(neq,maxk+2);
dif(:,1) = h * yp;

hinvGak = h * invGa(k);
nconhk = 0;                             % steps taken with current h and k

Miter = Mt - hinvGak * dfdy;

% Account for strongly state-dependent mass matrix.
if Mtype == 4
  psi = dif(:,K) * (G(K) * invGa(k));
  [dMpsidy,dMoptions.fac] = odenumjac(@odemxv, {Mfun,t,y,psi,Margs{:}}, Mt*psi, ...    
                                      dMoptions); %#ok<CCAT>
  Miter = Miter + dMpsidy;
end

% Use explicit scaling of the equations when solving DAEs.
if DAE
  RowScale = 1 ./ max(abs(Miter),[],2);
  Miter = sparse(one2neq,one2neq,RowScale) * Miter;
end
if issparse(Miter)
  [L,U,P,Q,R] = lu(Miter);
else
  [L,U,p] = lu(Miter,'vector');  
end  
ndecomps = ndecomps + 1;                
havrate = false;

% Allocate memory if we're generating output.
nout = 0;
tout = []; yout = [];
if nargout > 0
  if output_sol
    chunk = min(max(100,50*refine), refine+floor((2^11)/neq));      
    tout = zeros(1,chunk);
    yout = zeros(neq,chunk);
    kvec = zeros(1,chunk);
    dif3d = zeros(neq,maxk+2,chunk);
  else      
    if ntspan > 2                         % output only at tspan points
      tout = zeros(1,ntspan);
      yout = zeros(neq,ntspan);
    else                                  % alloc in chunks
      chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
      tout = zeros(1,chunk);
      yout = zeros(neq,chunk);
    end
  end  
  nout = 1;
  tout(nout) = t;
  yout(:,nout) = y;  
end

% Initialize the output function.
if haveOutputFcn
  feval(outputFcn,[t tfinal],y(outputs),'init',outputArgs{:});
end
  
% Cleanup the main ode function call
if ~isempty(odeArgs)
    ode = @(t,y) ode(t,y,odeArgs{:});
end

% THE MAIN LOOP

done = false;
at_hmin = false;
while ~done
  
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));
  if absh == hmin
    if at_hmin
      absh = abshlast;  % required by stepsize recovery
    end  
    at_hmin = true;
  else
    at_hmin = false;
  end  
  h = tdir * absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  if (absh ~= abshlast) || (k ~= klast)
    difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
    dif(:,K) = dif(:,K) * difRU(K,K);

    hinvGak = h * invGa(k);
    nconhk = 0;
    Miter = Mt - hinvGak * dfdy;
    if Mtype == 4
      Miter = Miter + dMpsidy;
    end    
    if DAE
      RowScale = 1 ./ max(abs(Miter),[],2);
      Miter = sparse(one2neq,one2neq,RowScale) * Miter;
    end
    if issparse(Miter)
      [L,U,P,Q,R] = lu(Miter);
    else  
      [L,U,p] = lu(Miter,'vector');
    end  
    ndecomps = ndecomps + 1;            
    havrate = false;
  end

  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true                            % Evaluate the formula.
    
    gotynew = false;                    % is ynew evaluated yet?
    while ~gotynew

      % Compute the constant terms in the equation for ynew.
      psi = dif(:,K) * (G(K) * invGa(k));

      % Predict a solution at t+h.
      tnew = t + h;
      if done
        tnew = tfinal;   % Hit end point exactly.
      end
      h = tnew - t;      % Purify h.
      pred = y + sum(dif(:,K),2);
      ynew = pred;
      
      % The difference, difkp1, between pred and the final accepted 
      % ynew is equal to the backward difference of ynew of order
      % k+1. Initialize to zero for the iteration to compute ynew.
      difkp1 = zeros(neq,1); 
      if normcontrol
        normynew = norm(ynew);
        invwt = 1 / max(max(normy,normynew),threshold);
        minnrm = 100*eps*(normynew * invwt);
      else
        invwt = 1 ./ max(max(abs(y),abs(ynew)),threshold);
        minnrm = 100*eps*norm(ynew .* invwt,inf);
      end

      % Mtnew is required in the RHS function evaluation.
      if Mtype == 2  % state-independent
        if odeIsFuncHandle
          Mtnew = Mfun(tnew,Margs{:}); % mass(t,p1,p2...)
        else                                     
          Mtnew = Mfun(tnew,ynew,Margs{:}); % mass(t,y,'mass',p1,p2...)
        end
      end
      
      % Iterate with simplified Newton method.
      tooslow = false;
      for iter = 1:maxit
        if Mtype >= 3 
          Mtnew = Mfun(tnew,ynew,Margs{:}); % state-dependent
        end
        rhs = hinvGak*ode(tnew,ynew) -  Mtnew*(psi+difkp1);
        if DAE                          % Account for row scaling.
          rhs = RowScale .* rhs;
        end

        % use overloaded subfunction mldivide which throws no warning.
        if issparse(Miter)
          del = Q * (U \ (L \ (P * (R \ rhs))));
        else  
          del = U \ (L \ rhs(p));
        end  

        if normcontrol
          newnrm = norm(del) * invwt;
        else
          newnrm = norm(del .* invwt,inf);
        end
        difkp1 = difkp1 + del;
        ynew = pred + difkp1;
        
        if newnrm <= minnrm
          gotynew = true;
          break;
        elseif iter == 1
          if havrate
            errit = newnrm * rate / (1 - rate);
            if errit <= 0.05*rtol       % More stringent when using old rate.
              gotynew = true;
              break;
            end
          else
            rate = 0;
          end
        elseif newnrm > 0.9*oldnrm
          tooslow = true;
          break;
        else
          rate = max(0.9*rate, newnrm / oldnrm);
          havrate = true;                 
          errit = newnrm * rate / (1 - rate);
          if errit <= 0.5*rtol             
            gotynew = true;
            break;
          elseif iter == maxit            
            tooslow = true;
            break;
          elseif 0.5*rtol < errit*rate^(maxit-iter)
            tooslow = true;
            break;
          end
        end
        
        oldnrm = newnrm;
      end                               % end of Newton loop
      nfevals = nfevals + iter;         
      nsolves = nsolves + iter;         
      
      if tooslow
        nfailed = nfailed + 1;          
        % Speed up the iteration by forming new linearization or reducing h.
        if ~Jcurrent || ~Mcurrent
          if ~Jcurrent  
            if Janalytic
              dfdy = Jac(t,y,Jargs{:});
            else
              f0 = ode(t,y);
              [dfdy,Joptions.fac,nF] = odenumjac(ode, {t,y}, f0, Joptions);
              nfevals = nfevals + nF + 1; 
            end             
            npds = npds + 1;            
            Jcurrent = true;
          end
          if ~Mcurrent
            Mt = Mfun(t,y,Margs{:});
            Mcurrent = true;
            if Mtype == 4
              [dMpsidy,dMoptions.fac] = odenumjac(@odemxv, {Mfun,t,y,psi,Margs{:}}, Mt*psi, ...
                                                  dMoptions); %#ok<CCAT>
            end
          end                       
        elseif absh <= hmin
          warning(message('MATLAB:ode15s:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));         
          solver_output = odefinalize(solver_name, sol,...
                                      outputFcn, outputArgs,...
                                      printstats, [nsteps, nfailed, nfevals,...
                                                   npds, ndecomps, nsolves],...
                                      nout, tout, yout,...
                                      haveEventFcn, teout, yeout, ieout,...
                                      {kvec,dif3d,idxNonNegative});
          if nargout > 0
            varargout = solver_output;
          end  
          return;
        else
          abshlast = absh;
          absh = max(0.3 * absh, hmin);
          h = tdir * absh;
          done = false;

          difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
          dif(:,K) = dif(:,K) * difRU(K,K);
          
          hinvGak = h * invGa(k);
          nconhk = 0;
        end
        Miter = Mt - hinvGak * dfdy;
        if Mtype == 4
          Miter = Miter + dMpsidy;
        end
        if DAE
          RowScale = 1 ./ max(abs(Miter),[],2);
          Miter = sparse(one2neq,one2neq,RowScale) * Miter;
        end
        if issparse(Miter)
          [L,U,P,Q,R] = lu(Miter);
        else  
          [L,U,p] = lu(Miter,'vector');
        end  
        ndecomps = ndecomps + 1;        
        havrate = false;
      end   
    end     % end of while loop for getting ynew
    
    % difkp1 is now the backward difference of ynew of order k+1.
    if normcontrol
      err = (norm(difkp1) * invwt) * erconst(k);
    else
      err = norm(difkp1 .* invwt,inf) * erconst(k);
    end
    if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
      if normcontrol
        errNN = norm( max(0,-ynew(idxNonNegative)) ) * invwt;
      else
        errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);
      end
      if errNN > rtol
        err = errNN;
      end
    end
    
    if err > rtol                       % Failed step
      nfailed = nfailed + 1;            
      if absh <= hmin
        warning(message('MATLAB:ode15s:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
        solver_output = odefinalize(solver_name, sol,...
                                    outputFcn, outputArgs,...
                                    printstats, [nsteps, nfailed, nfevals,...
                                                 npds, ndecomps, nsolves],...
                                    nout, tout, yout,...
                                    haveEventFcn, teout, yeout, ieout,...
                                    {kvec,dif3d,idxNonNegative});          
        if nargout > 0
          varargout = solver_output;
        end  
        return;
      end
      
      abshlast = absh;
      if nofailed
        nofailed = false;
        hopt = absh * max(0.1, 0.833*(rtol/err)^(1/(k+1))); % 1/1.2
        if k > 1
          if normcontrol
            errkm1 = (norm(dif(:,k) + difkp1) * invwt) * erconst(k-1);
          else
            errkm1 = norm((dif(:,k) + difkp1) .* invwt,inf) * erconst(k-1);
          end
          hkm1 = absh * max(0.1, 0.769*(rtol/errkm1)^(1/k)); % 1/1.3
          if hkm1 > hopt
            hopt = min(absh,hkm1);      % don't allow step size increase
            k = k - 1;
            K = 1:k;
          end
        end
        absh = max(hmin, hopt);
      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      if absh < abshlast
        done = false;
      end
      
      difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
      dif(:,K) = dif(:,K) * difRU(K,K);
      
      hinvGak = h * invGa(k);
      nconhk = 0;
      Miter = Mt - hinvGak * dfdy;
      if Mtype == 4
        Miter = Miter + dMpsidy;
      end      
      if DAE
        RowScale = 1 ./ max(abs(Miter),[],2);
        Miter = sparse(one2neq,one2neq,RowScale) * Miter;
      end
      if issparse(Miter)
        [L,U,P,Q,R] = lu(Miter);
      else   
        [L,U,p] = lu(Miter,'vector');
      end
      ndecomps = ndecomps + 1;          
      havrate = false;
      
    else                                % Successful step
      break;
      
    end
  end % while true
  nsteps = nsteps + 1;                  
  
  dif(:,k+2) = difkp1 - dif(:,k+1);
  dif(:,k+1) = difkp1;
  for j = k:-1:1
    dif(:,j) = dif(:,j) + dif(:,j+1);
  end
  
  NNreset_dif = false;
  if nonNegative && any(ynew(idxNonNegative) < 0)
    NNidx = idxNonNegative(ynew(idxNonNegative) < 0); % logical indexing
    ynew(NNidx) = 0;
    if normcontrol
      normynew = norm(ynew);
    end
    NNreset_dif = true;
  end   
  
  if haveEventFcn
    [te,ye,ie,valt,stop] = odezero(@ntrp15s,eventFcn,eventArgs,valt,...
                                   t,y,tnew,ynew,t0,h,dif,k,idxNonNegative);
    if ~isempty(te)
      if output_sol || (nargout > 2)
        teout = [teout, te]; %#ok<AGROW>
        yeout = [yeout, ye]; %#ok<AGROW>
        ieout = [ieout, ie]; %#ok<AGROW>
      end
      if stop               % Stop on a terminal event. 
        % Adjust the interpolation data to [t te(end)].                 
        taux = te(end) - (0:k)*(te(end) - t);
        yaux = ntrp15s(taux,t,y,tnew,ynew,h,dif,k,idxNonNegative);
        for j=2:k+1
          yaux(:,j:k+1) = yaux(:,j-1:k) - yaux(:,j:k+1);
        end
        dif(:,1:k) = yaux(:,2:k+1);        
        tnew = te(end);
        ynew = ye(:,end);
        h = tnew - t;
        done = true;
      end
    end
  end

  if output_sol
    nout = nout + 1;
    if nout > length(tout)
      tout = [tout, zeros(1,chunk)]; %#ok<AGROW>  requires chunk >= refine
      yout = [yout, zeros(neq,chunk)]; %#ok<AGROW> 
      kvec = [kvec, zeros(1,chunk)]; %#ok<AGROW>
      dif3d = cat(3,dif3d, zeros(neq,maxk+2,chunk));
    end
    tout(nout) = tnew; %#ok<AGROW>
    yout(:,nout) = ynew; %#ok<AGROW>
    kvec(nout) = k; %#ok<AGROW>
    dif3d(:,:,nout) = dif; %#ok<AGROW>
  end   
  
  if output_ty || haveOutputFcn 
    switch outputAt
     case 'SolverSteps'        % computed points, no refinement
      nout_new = 1;
      tout_new = tnew;
      yout_new = ynew;
     case 'RefinedSteps'       % computed points, with refinement
      tref = t + (tnew-t)*S;
      nout_new = refine;
      tout_new = [tref, tnew];
      yout_new = [ntrp15s(tref,[],[],tnew,ynew,h,dif,k,idxNonNegative), ynew];
     case 'RequestedPoints'    % output only at tspan points
      nout_new =  0;
      tout_new = [];
      yout_new = [];
      while next <= ntspan  
        if tdir * (tnew - tspan(next)) < 0
          if haveEventFcn && stop     % output tstop,ystop
            nout_new = nout_new + 1;
            tout_new = [tout_new, tnew]; %#ok<AGROW>
            yout_new = [yout_new, ynew]; %#ok<AGROW>
          end
          break;
        end
        nout_new = nout_new + 1;       
        tout_new = [tout_new, tspan(next)]; %#ok<AGROW>
        if tspan(next) == tnew
          yout_new = [yout_new, ynew]; %#ok<AGROW>
        else  
          yout_new = [yout_new, ntrp15s(tspan(next),[],[],tnew,ynew,h,dif,k,...
              idxNonNegative)]; %#ok<AGROW>
        end  
        next = next + 1;
      end
    end
    
    if nout_new > 0
      if output_ty
        oldnout = nout;
        nout = nout + nout_new;
        if nout > length(tout)
          tout = [tout, zeros(1,chunk)]; %#ok<AGROW> requires chunk >= refine
          yout = [yout, zeros(neq,chunk)]; %#ok<AGROW>
        end
        idx = oldnout+1:nout;        
        tout(idx) = tout_new; %#ok<AGROW>
        yout(:,idx) = yout_new; %#ok<AGROW>
      end
      if haveOutputFcn
        stop = feval(outputFcn,tout_new,yout_new(outputs,:),'',outputArgs{:});
        if stop
          done = true;
        end  
      end     
    end  
  end
  
  if done
    break
  end

  klast = k;
  abshlast = absh;
  nconhk = min(nconhk+1,maxk+2);
  if nconhk >= k + 2
    temp = 1.2*(err/rtol)^(1/(k+1));
    if temp > 0.1
      hopt = absh / temp;
    else
      hopt = 10*absh;
    end
    kopt = k;
    if k > 1
      if normcontrol
        errkm1 = (norm(dif(:,k)) * invwt) * erconst(k-1);
      else
        errkm1 = norm(dif(:,k) .* invwt,inf) * erconst(k-1);
      end
      temp = 1.3*(errkm1/rtol)^(1/k);
      if temp > 0.1
        hkm1 = absh / temp;
      else
        hkm1 = 10*absh;
      end
      if hkm1 > hopt 
        hopt = hkm1;
        kopt = k - 1;
      end
    end
    if k < maxk
      if normcontrol
        errkp1 = (norm(dif(:,k+2)) * invwt) * erconst(k+1);
      else
        errkp1 = norm(dif(:,k+2) .* invwt,inf) * erconst(k+1);
      end
      temp = 1.4*(errkp1/rtol)^(1/(k+2));
      if temp > 0.1
        hkp1 = absh / temp;
      else
        hkp1 = 10*absh;
      end
      if hkp1 > hopt 
        hopt = hkp1;
        kopt = k + 1;
      end
    end
    if hopt > absh
      absh = hopt;
      if k ~= kopt
        k = kopt;
        K = 1:k;
      end
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  if NNreset_dif  
    % Used dif for unperturbed solution to select order and interpolate.  
    % In perturbing ynew, defined NNidx.  Use now to reset dif to move along 
    % constraint.
    dif(NNidx,:) = 0;      
  end
  if normcontrol
    normy = normynew;
  end
  Jcurrent = Jconstant;
  switch Mtype
  case {0,1}
    Mcurrent = true;                    % Constant mass matrix I or M.
  case 2
    % M(t) has already been evaluated at tnew in Mtnew.
    Mt = Mtnew;
    Mcurrent = true;
  case {3,4}  % state dependent
    % M(t,y) has not yet been evaluated at the accepted ynew.
    Mcurrent = false;
  end
  
end % while ~done

solver_output = odefinalize(solver_name, sol,...
                            outputFcn, outputArgs,...
                            printstats, [nsteps, nfailed, nfevals,...
                                         npds, ndecomps, nsolves],...
                            nout, tout, yout,...
                            haveEventFcn, teout, yeout, ieout,...
                            {kvec,dif3d,idxNonNegative});
if nargout > 0
  varargout = solver_output;
end
