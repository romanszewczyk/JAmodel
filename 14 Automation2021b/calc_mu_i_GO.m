%
% The MIT License (MIT)
%
% Copyright (c) 2020 Roman Szewczyk
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
% Script calculates the parameters a and Kan for materials with GO anisotropy
% Required information about easy and hard axis initial relative permeability and saturation magnetization
%
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C. , Atherton D. "Theory of ferromagnetic hysteresis" Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Ramesh A., Jiles D., Roderik J. "A model of anisotropic anhysteretic magnetization" IEEE Transactions on Magnetics 32 (1996) 4234–4236. [Google Scholar] [CrossRef]
% [3] Ramesh A., Jiles D., Bi Y. "Generalization of hysteresis modeling to anisotropic materials" Journal of Applied Physics  81 (1997) 5585–5587.
% [4] Szewczyk R. "Validation of the Anhysteretic Magnetization Model for Soft Magnetic Materials with Perpendicular Anisotropy" Materials  7 (2014) 5109-5116.
% [5] Baghel, A. P. S., Kulkarni, S. V. "Hysteresis modeling of the grain-oriented laminations with inclusion of crystalline and textured structure in a modified 
%                                        Jiles-Atherton model" Journal of Applied Physics, 113(4) (2013) 043908.
% [6] Szewczyk R. " "

clear all
clc

fprintf('Calculation of parameters a and Kan of anhysteretic curve\n');
fprintf('of magnetic material with GO anisotropy\n');
fprintf('on the base of parallel and perpendicular relative inital permeability\n');
fprintf('and saturation magnetization Ms\n\n');

mi0 = 4.*pi*1e-7;

% -------- Input data --------

Ms = 1./mi0;
mi0_0 = 10000;
mi0_pi2 = 8000;

% ----------------------------
fprintf('Easy axis relative inital permeability = %1.4e\n', mi0_0);
fprintf('Hard axis relative inital permeability = %1.4e\n', mi0_pi2);
fprintf('Saturation magnetization Ms            = %1.4e (A/m)\n\n', Ms);

fprintf('Calculations...');

target= @(x) (mu_i_GO(Ms, x(1), x(2), 0) - mi0_0).^2 + ...
             (mu_i_GO(Ms, x(1), x(2), pi./2) - mi0_pi2).^2;
             % target function 

o = optimset("Display","notify");        
[res,fval] = fminsearch(target,[10 10],o);

fprintf(' Done.\n');

if (fval<1e-7) 
  fprintf('\nConvergence reached: fval= %1.2e \n\n', fval);
  fprintf('a=%2.4f (A/m), Kan= %2.4f (A/m) \n\n',res(1),res(2));
else
  fprintf('\n*** No convergence reached: fval= %1.2e *** \n\n', fval);
  fprintf('a=%2.4f (A/m), Kan= %2.4f (A/m) \n\n',res(1),res(2));
end
