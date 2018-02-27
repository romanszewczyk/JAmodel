function result = dMdH(a,k,c,Ms,alpha,M,H,Hstart,Hend,ModelType,IsoAniso,AnisoType,Kan,psi,IntType)
%
% The MIT License (MIT)
%
% Copyright (c) 2017 Roman Szewczyk, Peng Cheng
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
% dMdH is the function for solving main ODE (ordinary differential equation) determining Jiles-Atherton model
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl, Peng Cheng, Peng.Cheng@miun.se
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Chwastek K., Szczyglowski J. "Identification of a hysteresis model parameters with genetic algorithms" Mathematics and Computers in Simulation 71 (2006) 206.
%
% USAGE:
% result = dMdH(a,k,c,Ms,alpha,Kan,psi,M,H,Hstart,Hend,FixedStep)
% 
% INPUT:
% Hstart, Hend - magnetizing field, A/m (start and stop, A/m)
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
%             3 - Cheng Peng proposition 
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
% result - value of differential function in Jiles-Atherton model (scalar, vector or matrix, size of H)
%

if (max(size(a))>1 || max(size(k)>1) || max(size(c)>1) || max(size(Ms)>1) ...
   || max(size(alpha)>1) || max(size(Hstart)>1) || max(size(Hend)>1) ...
   || max(size(Kan)>1) || max(size(psi)>1))
 
  fprintf('\n*** ERROR in dMdH: among a, k, c, Kan, psi, Ms, alpha, Hstart, Hend at least one is not a scalar. ***\n');
  result=zeros(size(H));
  return;
  
end

if (size(H,1)~=size(M,1)) || (size(H,2)~=size(M,2))
 
  fprintf('\n*** ERROR in dMdH: size(H) is not equal size(M) ***\n');
  result=zeros(size(H));
  return;
  
end

Malpha=Mah_aniso(a,Ms,Kan,psi,H+alpha.*M,AnisoType,IntType);    % calculates anhysteretic anisotropic magnetization
                                    % to speed up further calculations

dM1 = (Malpha-M);                   % first part of equation


% Different types of Jiles-Atherton Model

if ModelType==0   % dMah/dHe

   if Hend>Hstart
      dM1(dM1<0)=0;
   end

   if Hend<Hstart
      dM1(dM1>0)=0;
   end                                 % check physical condition as it is explained in publication [2]

   if Hend<Hstart
      delta=-1;
   else
      delta=1;
   end         % end of if             % calculate delta - increase or decrease H


   dM2 = (1+c).*(delta.*k-alpha.*(Malpha-M));  % second part of equation

   dM3 = c./(1+c).*dMah_aniso(a,Ms,Kan,psi,H+alpha.*M,AnisoType,IntType);  % third part of equation

   dM2(abs(dM2)<1e-70)=1e-70;

   result=dM1./dM2+dM3;                        % calculate the equation [1]


end % of ModelType==0
   
   
if ModelType==1   % dMah/dH
   dM1_=dM1;
   
   if Hend>Hstart
      dM1(dM1<0)=0;
   end

   if Hend<Hstart
      dM1(dM1>0)=0;
   end                                 % check physical condition as it is explained in publication [2]

   if Hend<Hstart
      delta=-1;
   else
      delta=1;
   end         % end of if             % calculate delta - increase or decrease H

   dM2 = dM1./(delta.*k-alpha.*dM1_)+c.*dMah_aniso(a,Ms,Kan,psi,H+alpha.*M,AnisoType,IntType);  % second part of equation

   dM3 = 1+c-c.*alpha.*dMah_aniso(a,Ms,Kan,psi,H+alpha.*M,AnisoType,IntType);  % third part of equation

   dM3(abs(dM3)<1e-70)=1e-70;
   
   result=dM2./dM3;  % calculate the equation [1]

end % of ModelType==1  
   
   
if ModelType==2    % Venkataraman

   if Hend>Hstart
      dM1(dM1<0)=0;
   end

   if Hend<Hstart
      dM1(dM1>0)=0;
   end                                 % check physical condition as it is explained in publication [2]

   if Hend<Hstart
      delta=-1;
   else
      delta=1;
   end         % end of if             % calculate delta - increase or decrease H


   dM2 = delta.*k.*c.*dMah_aniso(a,Ms,Kan,psi,H+alpha.*M,AnisoType,IntType)+dM1;  % second part of equation

   dM3 = delta.*k-alpha.*dM1-delta.*k.*c.*alpha.*dMah_aniso(a,Ms,Kan,psi,H+alpha.*M,AnisoType,IntType);  % third part of equation
   
   dM3(abs(dM3)<1e-70)=1e-70;
   result=dM2./dM3;                        % calculate the equation [1]
  

end % of ModelType==2
 
   
if ModelType==3    % Cheng Peng

   if Hend>Hstart
      dM1(dM1<0)=0;
   end

   if Hend<Hstart
      dM1(dM1>0)=0;
   end                                 % check physical condition as it is explained in publication [2]

   if Hend<Hstart
      delta=-1;
   else
      delta=1;
   end         % end of if             % calculate delta - increase or decrease H


   dM2 = dM1.*(1-c)+c.*delta.*k.*dMah_aniso(a,Ms,Kan,psi,H+alpha.*M,AnisoType,IntType)+dM1;  % second part of equation

   dM3 = delta.*k-alpha.*(1-c).*(Malpha-M);  % third part of equation
   
   dM3(abs(dM3)<1e-70)=1e-70;
   
   result=dM2./dM3;                        % calculate the equation [1]
  

end % of ModelType==3



end  % of function

