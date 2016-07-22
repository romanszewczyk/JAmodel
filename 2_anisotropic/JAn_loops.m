function BsimT = JAn_loops(a,k,c,Ms,alpha,Kan,psi,HmeasT,SolverType,FixedStep)
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
% JAsingle_loop splits the given magnetizing field H into monotonous sections and control solving
% ODE stating Jiles-Atherton model
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Szewczyk R. "Computational problems connected with Jiles-Atherton model of magnetic hysteresis". Advances in Intelligent Systems and Computing (Springer) 267 (2014) 275.
% USAGE:
% [H,M] = JAsolver(a,k,c,Ms,alpha,Hstart,Hend,M0,SolverType)
% 
% IMPORTANT: H must be a COLUMN vector 
%
% INPUT:
% a    - quantifies domain density, A/m (scalar)
% k    - quantifies average energy required to break pinning side, A/m (scalar)
% c    - magnetization reversibility, 0..1 (scalar)
% Ms   - saturation magnetization, A/m (scalar)
% alpha - Bloch coefficient (scalar)
% Kan  - average magnetic anisotropy density, K/m3 (scalar)
% psi  - the angle between a direction of the easy axis of magnetic anisotropy and the direction of the magnetizing field, rad (scalar)
% HmeasT - magnetizing field, A/m (set of column vectors - matrix, A/m)
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
% Hw - set of output values of magnetizing field H, A/m (vector)
% Bw - set of output values of flux density B, T (vector)

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
 
 [Hw,Bw] = JAsingle_loop(a,k,c,Ms,alpha,Kan,psi,H,M0,SolverType,FixedStep);
  
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