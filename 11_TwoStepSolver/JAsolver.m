function [H,M] = JAsolver(a,k,c,Ms,alpha,Hstart,Hend,M0,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType)
% The MIT License (MIT)
%
% Copyright (c) 2016 Roman Szewczyk
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
% DESCRIPTION:
% JAsolver solves main ODE (ordinary differential equation) determining Jiles-Atherton model using different methods
% Important: JAsolver requires monotonous changes in the magnetizing field H
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Szewczyk R. "Computational problems connected with Jiles-Atherton model of magnetic hysteresis". Advances in Intelligent Systems and Computing (Springer) 267 (2014) 275.
%
% USAGE:
% [H,M] = JAsolver(a,k,c,Ms,alpha,Hstart,Hend,M0,SolverType)
% 
% IMPORTANT: H must be a scalar or COLUMN vector 
%
% INPUT:
% a    - quantifies domain density, A/m (scalar)
% k    - quantifies average energy required to break pinning side, A/m (scalar)
% c    - magnetization reversibility, 0..1 (scalar)
% Ms   - saturation magnetization, A/m (scalar)
% alpha - Bloch coefficient (scalar)
% Kan  - average magnetic anisotropy density, K/m3 (scalar)
% psi  - the angle between a direction of the easy axis of magnetic anisotropy and the direction of the magnetizing field, rad (scalar)
% H    - magnetizing field, A/m (scalar or column vector)
% Hstart - minimal value of magnetization, A/m (scalar)
% Hend - maximal value of magnetization, A/m (scalar)
% M0   - initial value of magnetization for ODE, A/m (scalar)
% SolverType - select the solver for ODE
%               1 - ode23()
%               2 - ode45()
%               3 - ode23s()
%               4 - rk4() 
% FixedStep - select the solver for integration:
%               1: quadtrapz(), 0: quadgk()
%
% OUTPUT:
% H - set of output values of magnetizing field H, A/m (vector)
% M - set of output values of magnetization M, A/m (vector)


if Hend==Hstart
   H=[Hend Hstart]';
   M=[M0 M0]';
   return
   end

options=odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',abs(Hend-Hstart)./10,'InitialStep',(Hend-Hstart)./10);    
  
dMdH_=@(H,M) dMdH(a,k,c,Ms,alpha,M,H,Hstart,Hend,ModelType,IsoAniso,AnisoType,Kan,psi,IntType);

switch (SolverType)

  case 1
  try
  [H,M] = ode23(dMdH_,[Hstart Hend],M0,options);
  catch
  H=[Hstart Hend];
  M=[0 0];
  fprintf('x');
  end_try_catch


  case 2
  [H,M] = ode45(dMdH_,[Hstart Hend],M0,options);
  
  case 3
  [H,M] = ode23s(dMdH_,[Hstart Hend],M0,options);
  
  case 4
  [H,M] = rk4(dMdH_,[Hstart Hend],M0,50);

  otherwise
  
  try
  [H,M] = ode23(dMdH_,[Hstart Hend],M0,options);
  catch
  H=[Hstart Hend];
  M=[0 0];
  fprintf('x');
  end_try_catch

endswitch
  
end
