function varargout = ode23(ode,tspan,y0,options,varargin)
%ODE23  Solve non-stiff differential equations, low order method.
%   [TOUT,YOUT] = ODE23(ODEFUN,TSPAN,Y0) integrates the system of
%   differential equations y' = f(t,y) from time TSPAN(1) to TSPAN(end)
%   with initial conditions Y0. Each row in the solution array YOUT
%   corresponds to a time in the column vector TOUT. 
%     * ODEFUN is a function handle. For a scalar T and a vector Y,
%       ODEFUN(T,Y) must return a column vector corresponding to f(t,y).
%     * TSPAN is a two-element vector [T0 TFINAL] or a vector with
%       several time points [T0 T1 ... TFINAL]. If you specify more than
%       two time points, ODE23 returns interpolated solutions at the
%       requested times.
%     * YO is a column vector of initial conditions, one for each equation.
%
%   [TOUT,YOUT] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS) specifies integration
%   option values in the fields of a structure, OPTIONS. Create the
%   options structure with <a href="matlab:helpview('matlab','MATLAB_HELP_ODESET')">odeset</a>.
%
%   [TOUT,YOUT,TE,YE,IE] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS) produces
%   additional outputs for events. An event occurs when a specified function
%   of T and Y is equal to zero. See <a href="matlab:helpview('matlab','MATLAB_HELP_EVENTS')">ODE Event Location</a> for details.
%
%   SOL = ODE23(...) returns a solution structure instead of numeric
%   vectors. Use SOL as an input to DEVAL to evaluate the solution at
%   specific points. Use it as an input to ODEXTEND to extend the
%   integration interval.
%
%   ODE23 can solve problems M(t,y)*y' = f(t,y) with mass matrix M that is
%   nonsingular. Use ODESET to set the 'Mass' property to a function handle
%   or the value of the mass matrix. ODE15S and ODE23T can solve problems
%   with singular mass matrices.
%
%   ODE23, ODE45, ODE78, and ODE89 are all single-step solvers that use
%   explicit Runge-Kutta formulas of different orders to estimate the error
%   in each step.
%     * ODE45 is for general use.
%     * ODE23 is useful for moderately stiff problems.
%     * ODE78 and ODE89 may be more efficient than ODE45 on non-stiff problems
%       that are smooth except possibly for a few isolated discontinuities.
%     * ODE89 may be more efficient than ODE78 on very smooth problems, when 
%       integrating over long time intervals, or when tolerances are tight.
%
%   Example
%         [t,y]=ode23(@vdp1,[0 20],[2 0]);   
%         plot(t,y(:,1));
%     solves the system y' = vdp1(t,y), using the default relative error
%     tolerance 1e-3 and the default absolute tolerance of 1e-6 for each
%     component, and plots the first component of the solution. 
%
%   Class support for inputs TSPAN, Y0, and the result of ODEFUN(T,Y):
%     float: double, single
%
%   See also ODE45, ODE78, ODE89, ODE113, ODE15S, ODE23S, ODE23T, ODE23TB,
%            ODE15I, ODESET, ODEPLOT, ODEPHAS2, ODEPHAS3, ODEPRINT, DEVAL,
%            ODEEXAMPLES, FUNCTION_HANDLE.

%   ODE23 is an implementation of the explicit Runge-Kutta (2,3) pair of
%   Bogacki and Shampine called BS23.  It uses a "free" interpolant of order
%   3.  Local extrapolation is done.

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-14-94
%   Copyright 1984-2021 The MathWorks, Inc.

solver_name = 'ode23';

% Check inputs
if nargin < 4
  options = [];
  if nargin < 3
    y0 = [];
    if nargin < 2
      tspan = [];
      if nargin < 1
        error(message('MATLAB:ode23:NotEnoughInputs'));
      end  
    end
  end
end

% Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0; 

[ode, odeIsFuncHandle, odeTreatAsMFile] = packageAsFuncHandle(ode);

% Output
output_sol = (odeIsFuncHandle && (nargout==1));      % sol = odeXX(...)
output_ty  = (~output_sol && (nargout > 0));  % [t,y,...] = odeXX(...)
% There might be no output requested...

sol = []; f3d = []; 
if output_sol
  sol.solver = solver_name;
  sol.extdata.odefun = ode;
  sol.extdata.options = options;                       
  sol.extdata.varargin = varargin;  
end  

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, ...
 options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType] = ...
    odearguments(odeIsFuncHandle, odeTreatAsMFile, solver_name, ode, tspan, y0, options, varargin);
nfevals = nfevals + 1;

% Handle the output
if nargout > 0
  outputFcn = odeget(options,'OutputFcn',[],'fast');
else
  outputFcn = odeget(options,'OutputFcn',@odeplot,'fast');
end
outputArgs  = {};      
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
[Mtype, M, Mfun] =  odemass(odeIsFuncHandle,ode,t0,y0,options,varargin);
if Mtype > 0  % non-trivial mass matrix
  Msingular = odeget(options,'MassSingular','no','fast');
  if strcmp(Msingular,'maybe')    
    warning(message('MATLAB:ode23:MassSingularAssumedNo'));
  elseif strcmp(Msingular,'yes')
    error(message('MATLAB:ode23:MassSingularYes'));
  end
  % Incorporate the mass matrix into ode and odeArgs.
  [ode,odeArgs] = odemassexplicit(odeIsFuncHandle,Mtype,ode,odeArgs,Mfun,M);
  f0 = ode(t0,y0,odeArgs{:});
  nfevals = nfevals + 1;
end

% Non-negative solution components
idxNonNegative = odeget(options,'NonNegative',[],'fast');
nonNegative = ~isempty(idxNonNegative);
if nonNegative  % modify the derivative function
  [ode,thresholdNonNegative] = odenonnegative(ode,y0,threshold,idxNonNegative);
  f0 = ode(t0,y0,odeArgs{:});
  nfevals = nfevals + 1;
end

t = t0;
y = y0;

% Allocate memory if we're generating output.
nout = 0;
tout = []; yout = [];
if nargout > 0
  if output_sol
    chunk = min(max(100,50*refine), refine+floor((2^11)/neq));      
    tout = zeros(1,chunk,dataType);
    yout = zeros(neq,chunk,dataType);
    f3d  = zeros(neq,4,chunk,dataType);
  else      
    if ntspan > 2                         % output only at tspan points
      tout = zeros(1,ntspan,dataType);
      yout = zeros(neq,ntspan,dataType);
    else                                  % alloc in chunks
      chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
      tout = zeros(1,chunk,dataType);
      yout = zeros(neq,chunk,dataType);
    end
  end  
  nout = 1;
  tout(nout) = t;
  yout(:,nout) = y;  
end

% Initialize method parameters.
pow = 1/3;
A = [1/2, 3/4, 1];
B = [
    1/2         0               2/9
    0           3/4             1/3
    0           0               4/9
    0           0               0
    ]; 
E = [-5/72; 1/12; 1/9; -1/8];
f = zeros(neq,4,dataType);
hmin = 16*eps(t);
if isempty(htry)
  % Compute an initial step size h using y'(t).
  absh = min(hmax, htspan);
  if normcontrol
    rh = (norm(f0) / max(normy,threshold)) / (0.8 * rtol^pow);
  else
    rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
  end
  if absh * rh > 1
    absh = 1 / rh;
  end
  absh = max(absh, hmin);
else
  absh = min(hmax, max(hmin, htry));
end
f(:,1) = f0;

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
while ~done
  
  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = tdir * absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true
    hB = h * B;
    f(:,2) = ode(t+h*A(1),y+f*hB(:,1));
    f(:,3) = ode(t+h*A(2),y+f*hB(:,2));
    tnew = t + h*A(3);
    if done
      tnew = tfinal;   % Hit end point exactly.
    end
    h = tnew - t;      % Purify h.     
    
    ynew = y + f*hB(:,3);
    f(:,4) = ode(tnew,ynew);
    nfevals = nfevals + 3;              
    
    % Estimate the error.
    NNrejectStep = false;
    if normcontrol
      normynew = norm(ynew);
      errwt = max(max(normy,normynew),threshold);      
      err = absh * (norm(f * E) / errwt);
      if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
        errNN = norm( max(0,-ynew(idxNonNegative)) ) / errwt ;
        if errNN > rtol
          err = errNN;
          NNrejectStep = true;
        end        
      end            
    else
      err = absh * norm((f * E) ./ max(max(abs(y),abs(ynew)),threshold),inf);
      if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
        errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);
        if errNN > rtol
          err = errNN;
          NNrejectStep = true;
        end              
      end            
    end
    
    % Accept the solution only if the weighted error is no more than the
    % tolerance rtol.  Estimate an h that will yield an error of rtol on
    % the next step or the next try at taking this step, as the case may be,
    % and use 0.8 of this value to avoid failures.
    if err > rtol                       % Failed step
      nfailed = nfailed + 1;            
      if absh <= hmin
        warning(message('MATLAB:ode23:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));        
        solver_output = odefinalize(solver_name, sol,...
                                    outputFcn, outputArgs,...
                                    printstats, [nsteps, nfailed, nfevals],...
                                    nout, tout, yout,...
                                    haveEventFcn, teout, yeout, ieout,...
                                    {f3d,idxNonNegative});
        if nargout > 0
          varargout = solver_output;
        end  
        return;
      end
      
      if nofailed
        nofailed = false;
        if NNrejectStep
          absh = max(hmin, 0.5*absh);
        else
          absh = max(hmin, absh * max(0.5, 0.8*(rtol/err)^pow));
        end
      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      done = false;
      
    else                                % Successful step

      NNreset_f4 = false;
      if nonNegative && any(ynew(idxNonNegative)<0)
        ynew(idxNonNegative) = max(ynew(idxNonNegative),0);
        if normcontrol
          normynew = norm(ynew);
        end
        NNreset_f4 = true;
      end  
      
      break;
      
    end
  end
  nsteps = nsteps + 1;                  
  
  if haveEventFcn
    [te,ye,ie,valt,stop] = ...
        odezero(@ntrp23,eventFcn,eventArgs,valt,t,y,tnew,ynew,t0,h,f,idxNonNegative);
    if ~isempty(te)
      if output_sol || (nargout > 2)
        teout = [teout, te]; %#ok<AGROW>
        yeout = [yeout, ye]; %#ok<AGROW>
        ieout = [ieout, ie]; %#ok<AGROW>
      end
      if stop               % Stop on a terminal event.               
        % Adjust the interpolation data to [t te(end)].         
        
        % Update the derivatives using the interpolating polynomial.
        taux = t + (te(end) - t)*A;        
        [~,f(:,2:4)] = ntrp23(taux,t,y,[],[],h,f,idxNonNegative);        

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
      tout = [tout, zeros(1,chunk,dataType)]; %#ok<AGROW> requires chunk >= refine
      yout = [yout, zeros(neq,chunk,dataType)]; %#ok<AGROW>
      f3d  = cat(3,f3d,zeros(neq,4,chunk,dataType)); 
    end
    tout(nout) = tnew; %#ok<AGROW>
    yout(:,nout) = ynew; %#ok<AGROW>
    f3d(:,:,nout) = f; %#ok<AGROW>
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
      yout_new = [ntrp23(tref,t,y,[],[],h,f,idxNonNegative), ynew];
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
          yout_new = [yout_new, ntrp23(tspan(next),t,y,[],[],h,f,idxNonNegative)]; %#ok<AGROW>
        end
        next = next + 1;
      end
    end
    
    if nout_new > 0
      if output_ty
        oldnout = nout;
        nout = nout + nout_new;
        if nout > length(tout)
          tout = [tout, zeros(1,chunk,dataType)]; %#ok<AGROW> requires chunk >= refine
          yout = [yout, zeros(neq,chunk,dataType)]; %#ok<AGROW>
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
  
  % If there were no failures compute a new h.
  if nofailed
    % Note that absh may shrink by 0.8, and that err may be 0.
    temp = 1.25*(err/rtol)^pow;
    if temp > 0.2
      absh = absh / temp;
    else
      absh = 5.0*absh;
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  if normcontrol
    normy = normynew;
  end
  if NNreset_f4
    % Used f4 for unperturbed solution to interpolate.  
    % Now reset f4 to move along constraint.
    f(:,4) = ode(tnew,ynew);
    nfevals = nfevals + 1;
  end 
  f(:,1) = f(:,4);  % Already have f(tnew,ynew)
  
end

solver_output = odefinalize(solver_name, sol,...
                            outputFcn, outputArgs,...
                            printstats, [nsteps, nfailed, nfevals],...
                            nout, tout, yout,...
                            haveEventFcn, teout, yeout, ieout,...
                            {f3d,idxNonNegative});
if nargout > 0
  varargout = solver_output;
end  
