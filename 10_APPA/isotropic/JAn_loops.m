function BsimT = JAn_loops(HmeasT,a,k,c,Ms,alpha,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType)
%
% The MIT License (MIT)
%
% Copyright (c) 2017 Roman Szewczyk
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
% JAn_loop solves Jiles-Atherton model for n hysteresis loops 
%          (according to the numner of columns in HmeasT)
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1]  Roman Szewczyk, Peng Cheng "Open source implementation of different variants of
%      Jiles-Atherton model of magnetic hysteresis loops" Acta Physica Polonica A (submitted)
%
% USAGE:
% BsimT = JAn_loops(a,k,c,Ms,alpha,Hstart,Hend,M0,SolverType)
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
% SolverType - select the solver for ODE
%               0 - ode23() - DEFAULT
%               1 - ode45()
%               2 - ode23s()
%               3 - rk4() 
% 
% IntType - select the solver for integration:
%               0: quadtrapz() - DEFAULT
%               1: quadgk()
%
% OUTPUT:
% BsimT - set of output values of flux density B, T (vector)


% determine default values for not determined variables
if nargin<13 || isempty(IntType)
   IntType=0;
end
if nargin<12 || isempty(psi)
   psi=0;
end
if nargin<11 || isempty(Kan)
   Kan=0;
end
if nargin<10 || isempty(AnisoType)
   AnisoType=0;
end
if nargin<9 || isempty(IsoAniso)
   IsoAniso=0;
end
if nargin<8 || isempty(SolverType)
   SolverType=0;
end
if nargin<7 || isempty(ModelType)
   ModelType=0;
end

if ((a<0) || (k<0) || (c<0) || (c>1) || (alpha<0) || (Kan<0)) 
    fprintf('\n\n*** ERROR in JAn_loops: unphysical parameters.\n\n');
    BsimT=zeros(size(HmeasT));
    return
end           % check if parameters are physical


BsimT=zeros(size(HmeasT));                    % Prepare output variable

InitialCurve=0;
if sum(HmeasT(1,:))==0
    InitialCurve=1;
end

M0=0; % sample demagnetized at the beginning

for lT=1:size(HmeasT,2),

 H=HmeasT(:,lT);

 if InitialCurve==0
     H=[0; H];
 end                  % check if there is initial magnetization curve in H
 
 [Hw,Bw] = JAsingle_loop(H,a,k,c,Ms,alpha,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType);
  
 if InitialCurve==0
     Hw(1)=[];
     Bw(1)=[];
     H(1)=[];
     
     Bw=Bw+linspace(0,Bw(1)-Bw(numel(Bw)),numel(Bw))';       % Close hysteresis loop
     Bw=Bw-(max(Bw)+min(Bw))./2    ;                         % Symmetrize Bw
          
 end
 
 BsimT(:,lT)=Bw;          % Remember the result

end


end