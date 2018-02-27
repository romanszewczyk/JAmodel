function result = Mah_aniso(a, Ms, Kan, psi, HeT, AnisoType,IntType)
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
% Mah_aniso is the function for calculation of anhysteretic, ANISOTROPIC magnetization in Jiles-Atherton model
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C. , Atherton D. "Theory of ferromagnetic hysteresis" Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Ramesh A., Jiles D., Roderik J. "A model of anisotropic anhysteretic magnetization" IEEE Transactions on Magnetics 32 (1996) 4234–4236. [Google Scholar] [CrossRef]
% [3] Ramesh A., Jiles D., Bi Y. "Generalization of hysteresis modeling to anisotropic materials" Journal of Applied Physics  81 (1997) 5585–5587.
% [4] Szewczyk R. "Validation of the Anhysteretic Magnetization Model for Soft Magnetic Materials with Perpendicular Anisotropy" Materials  7 (2014) 5109-5116.
%
% USAGE:
% M = Mah_aniso(a, Ms, Kan, psi, He)
% 
% INPUT:
% a    - quantifies domain density, A/m (scalar)
% Ms   - saturation magnetization, A/m (scalar)
% Kan  - average magnetic anisotropy density, K/m3 (scalar)
% psi  - the angle between a direction of the easy axis of magnetic anisotropy and the direction of the magnetizing field, rad (scalar)
% HeT  - effective magnetizing field, A/m (scalar, vector. Matrix is not acceptable)
% FixedStep - alghoritm for numerical integration:
%	      0 - quadgk() - Gaus-Kronrod quadrature
% 	    1 - trapz()  - trapezoidal
%
% AnisoType - select the type of anisotropy model
%               0: uniaxial anisotropy
%               1: GO anisotropy
%
%
% OUTPUT:
% result - value of anhysteretic, anisotropic magnetization in Jiles-Atherton model, A/m (scalar, vector or matrix, size of He)
%

if (max(size(a))>1 || max(size(Ms)>1))
  
  fprintf('\n*** ERROR in Mah_aniso: a or Ms is not a scalar. ***\n');
  return;
  
end
  
if (min(size(HeT))>1)

  fprintf('\n*** ERROR in Mah_aniso: HeT is a matrix. ***\n');
  return;
  
end

if (Ms==0) || (a==0)

  result=zeros(size(He));
  return;

end  % of if

i=(HeT==0);        % find singularities for He==0

HeT(i)=1;          % remove singularities

if Kan==0
   result=Ms.*(coth(HeT./a)-(a./HeT));   % calculate anhysteretic, isotropic magnetization due to lack of anisotropy
else
   result=HeT.*0;

   for j=1:numel(HeT)

       He=HeT(j);

       if AnisoType==1  % GO anisotropy
       e1 = @(theta) (He)./a.*cos(theta)-Kan./(Ms.*a.*4.*pi.*1e-7).*((cos(psi-theta)).^2.*(sin(psi-theta)).^2+0.25.*(sin(psi-theta)).^4);
       e2 = @(theta) (He)./a.*cos(theta)-Kan./(Ms.*a.*4.*pi.*1e-7).*((cos(psi+theta)).^2.*(sin(psi+theta)).^2+0.25.*(sin(psi+theta)).^4);
       F1 = @(theta) exp(0.5.*(e1(theta)+e2(theta))).*sin(theta).*cos(theta);
       F2 = @(theta) exp(0.5.*(e1(theta)+e2(theta))).*sin(theta);
       else   % uniaxial anisotropy
       e1 = @(theta) (He)./a.*cos(theta)-Kan./(Ms.*a.*4.*pi.*1e-7).*(sin(psi-theta)).^2;
       e2 = @(theta) (He)./a.*cos(theta)-Kan./(Ms.*a.*4.*pi.*1e-7).*(sin(psi+theta)).^2;
       F1 = @(theta) exp(0.5.*(e1(theta)+e2(theta))).*sin(theta).*cos(theta);
       F2 = @(theta) exp(0.5.*(e1(theta)+e2(theta))).*sin(theta);       
       end
       
       
       if IntType==0
          p1=quadtrapz(F1,0,pi,90);
          p2=quadtrapz(F2,0,pi,90);
       else
          p1=quadgk(F1,0,pi);
          p2=quadgk(F2,0,pi);  
       end

       if ((p2==0) || isnan(p1) || isnan(p2))
          result(j) = 0; 			% something wrong happen
       else
          result(j) = Ms.*p1./p2; 
       end
    end % end of for


end   % of if


result(i)=0;      % correct values for singularities


end  % of function

