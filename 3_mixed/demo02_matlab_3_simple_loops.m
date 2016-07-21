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
% Demonstration of solving Jiles-Atherton model of three hysteresis loops 
% of mixed isotropic and anisotropic magnetic material, with increase of magnetizing field amplitude
% 
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Ramesh A., Jiles D., Roderik J. "A model of anisotropic anhysteretic magnetization" IEEE Transactions on Magnetics 32 (1996) 4234–4236. [Google Scholar] [CrossRef]
% [3] Ramesh A., Jiles D., Bi Y. "Generalization of hysteresis modeling to anisotropic materials" Journal of Applied Physics  81 (1997) 5585–5587.
% [4] Szewczyk R. "Validation of the Anhysteretic Magnetization Model for Soft Magnetic Materials with Perpendicular Anisotropy" Materials  7 (2014) 5109-5116.
%
% USAGE:
% demo02_matlab_3_simple_loops
% 


clear all
clc

fprintf('\n\nDemonstration of Jiles-Atherton model - three hysteresis loops.');
fprintf('\nDemonstration optimized for MATLAB. For OCTAVE please use demo02_octave_3_simple_loops.m ');
fprintf('\nDemonstration requires odepkg installed.\n\n');

% Calculate Jiles-Atheton model for the case presented in [1]
 
% prepare variables for modelling

Ms=1e6;
a=50;
alpha=1e-5;
k=30;
c=0.8;        
Kan=1e3;
psi=0;   
t=0.6;  % Parameters of Jiles-Atherton model specified in [1]

fprintf('Calculation for parameters: \na=%2.2f(A/m), k=%2.2f(A/m), c=%1.2f, Ms=%1.2e(A/m), alpha=%1.2e Kan=%1.2eJ/m3 psi=0 t=%2.2f \n\n',a,k,c,Ms,alpha,Kan,t);

H=[0:0.01:1 1:-0.01:-1 -1:0.01:1]'; % magnetizing field H - column vector
H=[H.*100 H.*200 H.*500];

SolverType=1; % ode23() solver
FixedStep=1;

fprintf('Calculations started. Expected time of calculations: below 3 minutes.\n');

fprintf('\nCalculations for parallel anisotropy...');
psi=0;
Bmodel_par = JAn_loops(a,k,c,Ms,alpha,Kan,psi,t,H,SolverType,FixedStep);
fprintf('done (red line).');

fprintf('\nCalculations for perpendicular anisotropy...');
psi=pi./2;
Bmodel_per = JAn_loops(a,k,c,Ms,alpha,Kan,psi,t,H,SolverType,FixedStep);
fprintf('done (blue line).');

fprintf('\nCalculations for isotropic material...');
Kan=0;
Bmodel_neutral = JAn_loops(a,k,c,Ms,alpha,Kan,psi,t,H,SolverType,FixedStep);
fprintf('done (black line).\n\n');

plot(H, Bmodel_par,'r',H, Bmodel_per,'b',H, Bmodel_neutral,'k');
xlabel('H (A/m)');
ylabel('B (T)');
grid;

