function MahT = Mah_loop(HmeasT,a,Ms,alpha,AnisoType,Kan,psi,IntType)
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
%Mah_loop calculates anhysteretic curve of Jiles-Atherton model for single hysteresis loop 
%         %
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1]  
%
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
if nargin<8 || isempty(IntType)
   IntType=0;
end
if nargin<7 || isempty(psi)
   psi=0;
end
if nargin<6 || isempty(Kan)
   Kan=0;
end
if nargin<5 || isempty(AnisoType)
   AnisoType=0;
end

if ((a<0) || (alpha<0) || (Kan<0)) 
    fprintf('\n\n*** ERROR in JAn_loops: unphysical parameters.\n\n');
    Mah_loop=zeros(size(HmeasT));
    return
end           % check if parameters are physical

if size(HmeasT,2)>1
    fprintf('\n\n*** ERROR HmeasT must be the column vector.\n\n');
    Mah_loop=zeros(size(HmeasT));
    return
end           % check if Hmeast is column vector



MahT=HmeasT.*0;                   % Prepare output variable



for lT=1:numel(HmeasT),

fun=@(M) M-Mah_aniso(a, Ms, Kan, psi, HmeasT(lT)+alpha.*M, AnisoType,IntType);

[x, fval, info, output] = fzero (fun,[Ms -1.*Ms]);

MahT(lT)=x;

end


end
