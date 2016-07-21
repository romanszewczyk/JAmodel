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
% Demonstration of solving simple case of Jiles-Atherton model of hysteresis loop
% of mixed isotropic and anisotropic magnetic material
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
% demo01_octave_simple_loop
% 
% IMPORTANT: Demo requires "odepkg" installed and loaded  
% "odepkg" is package for solving ordinary differential equations.
%

clear all
clc

page_screen_output(0);
page_output_immediately(1);  % print immediately at the screen


fprintf('\n\nSimple demonstration of Jiles-Atherton model for ANISOTROPIC materials.');
fprintf('\nDemonstration optimized for OCTAVE. For MATLAB please use demo01_matlab_simple_loop.m ');
fprintf('\nDemonstration requires odepkg installed.\n\n');


% check if odepkg is installed. Load odepkg if installed, but not loaded.

inst_pkg = pkg ("list");

[i,j]=size(inst_pkg);

odepkg_inst=0;
odepkg_loaded=0;

for i=1:j
    if size(findstr(inst_pkg{1,i}.name,'odepkg'))>0
       odepkg_inst=1;
       if inst_pkg{1,i}.loaded==1;
          odepkg_loaded=1;
       end
    end
end

if odepkg_inst==0
   fprintf('\n *** ERROR: odepkg must be installed to solve ODEs.\n To solve problem try: pkg install -forge odepkg\n\n');
   return
else
   fprintf('\n odepkg installed...ok.');
end
   
 if odepkg_loaded==0
   fprintf('\n WARNING: odepkg is installed but not loaded.\n');
   pkg load odepkg
   fprintf(' Problem solved: odepkg is loaded now.\n\n');
   else
   fprintf('\n odepkg loaded...ok.\n\n');
end


% Calculate Jiles-Atheton model for given case with parallel and perpendicular anisotropy, as well as for isotropic material
 
% prepare variables for modelling

Ms=1e6;
a=50;
alpha=1e-5;
k=50;
c=0.8;        
Kan=1e3;
psi=0;   
t=0.6 % Parameters of Jiles-Atherton model 

fprintf('Calculation for parameters: \na=%2.2f(A/m), k=%2.2f(A/m), c=%1.2f, Ms=%1.2e(A/m), alpha=%1.2e Kan=%1.2eJ/m3 psi=0 rad, t=%2.2f \n\n',a,k,c,Ms,alpha,Kan,t);

H=[0:50:500 500:-5:-500 -500:5:500]'; % magnetizing field H - column vector

M0=0;         % sample demagnetized at the beginning
SolverType=1; % ode23() solver
FixedStep=1;


fprintf('Calculations started. Expected time of calculations: below 1 minute.\n');

fprintf('\nCalculations for parallel anisotropy...');
psi=0;
[Hmodel_par,Bmodel_par] = JAsingle_loop(a,k,c,Ms,alpha,Kan,psi,t,H,M0,SolverType,FixedStep);
fprintf('done (red line).');

fprintf('\nCalculations for perpendicular anisotropy...');
psi=pi./2;
[Hmodel_per,Bmodel_per] = JAsingle_loop(a,k,c,Ms,alpha,Kan,psi,t,H,M0,SolverType,FixedStep);
fprintf('done (blue line).');

fprintf('\nCalculations for isotropic material...');
Kan=0;
[Hmodel_neutral,Bmodel_neutral] = JAsingle_loop(a,k,c,Ms,alpha,Kan,psi,t,H,M0,SolverType,FixedStep);
fprintf('done (black line).\n\n');

plot(Hmodel_par, Bmodel_par,'r',Hmodel_per, Bmodel_per,'b',Hmodel_neutral, Bmodel_neutral,'k');
xlabel('H (A/m)');
ylabel('B (T)');
grid;


