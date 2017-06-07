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
% Visualization of identification of parameters of Jiles-Atherton model of four hysteresis loops with increase of magnetizing field amplitude
% 
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] 
%
% USAGE:
% demo03_octave_simple_parametrs_identification
% 
% IMPORTANT: Demo requires "odepkg", "struct" and "optim" packages installed and loaded  
%

clear all
clc

page_screen_output(0);
page_output_immediately(1);  % print immediately at the screen


fprintf('\n\nDemonstration of identification of Jiles-Atherton models parameters for three hysteresis loops.');
fprintf('\nDemonstration optimized for OCTAVE. For MATLAB please use demo03_matlab_simple_parameters_identification.m ');
fprintf('\nDemonstration requires odepkg, struct and optim packages installed.\n\n');

% check if odepkg is installed. Load odepkg if installed, but not loaded.
ChkPkg('odepkg');

% check if struct is installed. Load odepkg if installed, but not loaded.
ChkPkg('struct');

% check if optim is installed. Load odepkg if installed, but not loaded.
ChkPkg('optim');

% Load measured B(H) characterisitcs of Mn-Zn ferrite

cd ('Characterisitcs_anisotropic_mat');
load('H_M130_27s.mat');
load('B_M130_27s.mat');
cd ('..');

HmeasT=HmeasT(:,3);
BmeasT=BmeasT(:,3);

fprintf('Load measured B(H) characterisitcs of M130-27s GO electrical steel... done\n\n');

% prepare starting point for optimisation

mi0=4.*pi.*1e-7;
psi=0;

load('demo03_results.mat');

%JApoint0=[1 1 1 1 1 1];
%JApoint_res=[20  12.2  0.79  1.1250e+006  7.5e-007  244];

ModelType=1;
SolverType=1;
IsoAniso=1;
AnisoType=1;
IntType=0;

func = @(JApointn) JAn_loops_target( JApointn, JApoint0, HmeasT, BmeasT, ModelType,SolverType,IsoAniso,AnisoType,IntType);

fprintf('\n\nCalculations...\n\n');

Ftarget=func(JApoint_res);

BsimT = JAn_loops(HmeasT,JApoint0(1).*JApoint_res(1), JApoint0(2).*JApoint_res(2), JApoint0(3).*JApoint_res(3), ...
 JApoint0(4).*JApoint_res(4), JApoint0(5).*JApoint_res(5), ModelType,SolverType,IsoAniso,AnisoType, JApoint0(6).*JApoint_res(6),psi,IntType);


fprintf('Results of optimisation:\n'); 
fprintf('Target function value: Ftarget=%f\n',Ftarget);
fprintf('JA model params: a=%f(A/m), k=%f(A/m), c=%f, Ms=%e(A/m), alpha=%e, Kan=%e (J/m3), psi=0\n',  ...
 JApoint0(1).*JApoint_res(1), JApoint0(2).*JApoint_res(2), JApoint0(3).*JApoint_res(3), ...
 JApoint0(4).*JApoint_res(4), JApoint0(5).*JApoint_res(5), JApoint0(6).*JApoint_res(6));
 
plot(HmeasT, BmeasT,'r',HmeasT,BsimT,'k');
xlabel('H (A/m)');
ylabel('B (T)');
grid;

r2=corr(BsimT,BmeasT);
r2=r2(1,1).^2;
fprintf('JA model params: R2=%f \n\n', r2);
