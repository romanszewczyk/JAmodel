function Ftarget = JAn_loops_target(JApointn,JApoint0,HmeasT,BmeasT,ModelType,SolverType,IsoAniso,AnisoType,IntType)

%
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
% JAn_loop_target_mod calculates target function for optimisation. 
% Target function is the sum of squared differences between results of modelling and results of measurements.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1]  Roman Szewczyk, Peng Cheng "Open source implementation of different variants of
%      Jiles-Atherton model of magnetic hysteresis loops" Acta Physica Polonica A (submitted)
%
% USAGE:
% [H,M] = JAn_loops(a,k,c,Ms,alpha,Hstart,Hend,M0,SolverType)
% 
% IMPORTANT: H must be a COLUMN vector 
%
% INPUT:
% HmeasT - magnetizing field, A/m (set of column vectors - matrix, A/m)
% a    - quantifies domain density, A/m (scalar)
% k    - quantifies average energy required to break pinning side, A/m (scalar)
% c    - magnetization reversibility, 0..1 (scalar)
% Ms   - saturation magnetization, A/m (scalar)
% alpha - Bloch coefficient (scalar)
%
% ModelType - variation of the Jiles-Atherton model (please refer to [1])
%             0 - dMah/dHe - DEFAULT
%             1 - dMah/dH 
%             2 - Venkataraman bulk magnetic hysteresis model 
%
% SolverType - select the solver for ODE
%               0 - ode23() - DEFAULT
%               1 - ode45()
%               2 - ode23s()
%               3 - rk4() 
%
% IsoAniso - isotropic or anisotropic model
%             0 - isotropic model - DEFAULT
%             1 - anisotropic model
%
% AnisoType - 
%
% Kan  - average magnetic anisotropy density, K/m3 (scalar)
% psi  - the angle between a direction of the easy axis of magnetic anisotropy and the direction of the magnetizing field, rad (scalar)
% 
% IntType - select the solver for integration:
%               0: quadtrapz() - DEFAULT
%               1: quadgk()
%
% OUTPUT:
% BsimT - set of output values of flux density B, T (vector)


if numel(JApointn)<6
  JApointn(6)=0;
  JApoint0(6)=0;
end
if numel(JApointn)<7
  JApointn(7)=0;
  JApoint0(7)=0;
end

% determine default values for not determined variables
if nargin<9 || isempty(IntType)
   IntType=0;
end
if nargin<8 || isempty(AnisoType)
   AnisoType=0;
end
if nargin<7 || isempty(IsoAniso)
   IsoAniso=0;
end
if nargin<6 || isempty(SolverType)
   SolverType=0;
end
if nargin<5 || isempty(ModelType)
   ModelType=0;
end

if size(HmeasT)~=size(BmeasT)
    fprintf('\n\n*** ERROR in JAn_loops_target: diferent sizes HmeasT and BmeasT.\n\n');
    Ftarget=1e30;
    return
end  

JApoint=JApointn.*JApoint0;

if numel(JApoint)~=7
    fprintf('\n\n*** ERROR in JAn_loops_target: wrong size of JApoint.\n\n');
    Ftarget=1e30;
    return
end  

a=JApoint(1);
k=JApoint(2);
c=JApoint(3);
Ms=JApoint(4);
alpha=JApoint(5);
Kan=JApoint(6);
psi=JApoint(7);

if ((a<0) || (k<0) || (c<0) || (c>1) || (alpha<0) || (Kan<0)) 
    fprintf('\n\n*** ERROR in JAn_loops_target: unphysical parameters.\n\n');
    Ftarget=1e30;
    return
end           % check if parameters are physical

BsimT = JAn_loops(HmeasT,a,k,c,Ms,alpha,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType);

Ftarget=sum(sum((BmeasT-BsimT).^2));  % calculate target function

if isnan(Ftarget)==1
    Ftarget=1e30;                   % check if the target function is NaN. In such a case make it very high
end 

end
