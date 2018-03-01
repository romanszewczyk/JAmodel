function result = dMah_aniso(a,Ms,Kan,psi,He,AnisoType,IntType)
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
% dMah_aniso is the function for numerical differentiation of anhysteretic, anisotropic magnetization over effective magnetic field He in Jiles-Atherton model
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% RELATED PUBLICATION(S):
% [1] Jiles D. C. , Atherton D. "Theory of ferromagnetic hysteresis" Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Ramesh A., Jiles D., Roderik J. "A model of anisotropic anhysteretic magnetization" IEEE Transactions on Magnetics 32 (1996) 4234–4236. [Google Scholar] [CrossRef]
% [3] Ramesh A., Jiles D., Bi Y. "Generalization of hysteresis modeling to anisotropic materials" Journal of Applied Physics  81 (1997) 5585–5587.
% [4] Szewczyk R. "Validation of the Anhysteretic Magnetization Model for Soft Magnetic Materials with Perpendicular Anisotropy" Materials  7 (2014) 5109-5116.
%
% USAGE:
% dM = dMah_aniso(a, Ms, Kan, psi, He, FixedStep)
% 
% INPUT:
% a    - quantifies domain density, A/m (scalar)
% Ms   - saturation magnetization, A/m (scalar)
% He   - effective magnetizing field, A/m (scalar, vector. Matrix is not acceptable.)
% Kan  - average magnetic anisotropy density, K/m3 (scalar)
% psi  - the angle between a direction of the easy axis of magnetic anisotropy and the direction of the magnetizing field, rad (scalar)
% HeT  - effective magnetizing field, A/m (scalar, vector. Matrix is not acceptable)
% FixedStep - alghoritm for numerical integration:
%	      0 - quadgk() - Gaus-Kronrod quadrature
% 	    1 - quadtrapz()  - trapezoidal
%
% OUTPUT:
% result - value of numerical differentiation of anhysteretic, anisotropic magnetization over effective magnetic field in Jiles-Atherton model, A/m (scalar, vector or matrix, size of He)
%

if (max(size(a))>1 || max(size(Ms)>1))
  
  fprintf('\n*** ERROR in dMah_aniso: a or Ms is not a scalar. ***\n');
  return

end

if min(size(He))>1

  fprintf('\n*** ERROR in dMah_aniso: HeT is a matrix. ***\n');
  return
  
end  
  

result=(Mah_aniso(a, Ms, Kan, psi, He+1e-6,AnisoType,IntType)-Mah_aniso(a, Ms, Kan, psi, He-1e-6,AnisoType,IntType))./2e-6;   % value of numerical differentiation

end  % of function

