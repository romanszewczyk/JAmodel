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
% Script calculates the parameter a for isotropic material
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C. , Atherton D. "Theory of ferromagnetic hysteresis" Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Szewczyk R. " "
%

clear all
clc

fprintf('Calculation of parameter a of anhysteretic curve\n');
fprintf('of isotropic magnetic material\n');
fprintf('on the base of its relative inital permeability\n');
fprintf('and saturation magnetization Ms\n\n');

mi0 = 4.*pi*1e-7;

% -------- Input data --------

Ms = 1.8./mi0;
mi_i = 1000;

% ----------------------------
fprintf('Relative inital permeability = %1.4e\n', mi_i);
fprintf('Saturation magnetization Ms  = %1.4e (A/m)\n\n', Ms);

fprintf('Calculations...');
     
res = Ms./(3.*mi_i);

fprintf(' Done.\n');

fprintf('\n');
fprintf('a=%2.4f (A/m) \n\n',res);
