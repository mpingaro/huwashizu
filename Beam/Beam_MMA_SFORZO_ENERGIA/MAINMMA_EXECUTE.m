%---------------------------------------------------------------------
%  This is the file beammain.m.  Version September 2007.
%  Written by Krister Svanberg <krille@math.kth.se>.
%
%  This file contains a main program for using MMA to solve
%  a problem defined by the users files beaminit.m
%  (which must be run before beammain.m) and beam2.m.
%
%%%% If outeriter=0, the user should now calculate function values
%%%% and gradients of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx:
close all
clc

addpath C:\Users\Paolo\Documents\Paolo\MMA\MATLAB

if outeriter < 0.5
  [f0val,df0dx,fval,dfdx] = ottigoal_3(xval);
  innerit=0;
  outvector1 = [outeriter innerit f0val fval']
  outvector2 = xval'
end
%
%%%% The iterations start:
kktnorm = kkttol+10;
f0val = f0valsoglia+50;
fbest = 1.e+9;
xbest = [];
strainbest = [];
stressbest = [];

outit = 0;
while kktnorm > kkttol & outit < maxoutit
%while f0val > f0valsoglia && outit < maxoutit && kktnorm > kkttol 
  outit   = outit+1;
  outeriter = outeriter+1;
%%%% The MMA subproblem is solved at the point xval:
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
  mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
  f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);

%%%% Some vectors are updated:
  xold2 = xold1;
  xold1 = xval;
  xval  = xmma;
%%%% The user should now calculate function values and gradients
%%%% of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx.
  [f0val,df0dx,fval,dfdx] = ottigoal_3(xval);
  if f0val<fbest
      fbest = f0val;
      xbest = xval;
      strainvest = strain;
      stressbest = stress;
  end
%%%% The residual vector of the KKT conditions is calculated:
  [residu,kktnorm,residumax] = ...
  kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
           xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
  outvector1 = [outeriter innerit f0val fval'];
  outvector2 = xval';
  [outit/maxoutit kktnorm/kkttol f0val norm(df0dx)]
%
end
save xval xval
 
%% PLOT SOLUTION

% Deformed mesh
nnod = size(coordinates,1);
Ux = spost(1:2:2*nnod); 
Uy = spost(2:2:2*nnod);
defo = coordinates + 0.1*[Ux, Uy];

plotsolution_goal(coordinates,element,xval);

plotsolution_goal(coordinates,element,xbest);

%plotsolution_andmises(coordinates,element,0.1*defo,strain,stress,xval);

%---------------------------------------------------------------------



