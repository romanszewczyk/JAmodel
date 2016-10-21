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
% View of identified parameters of Jiles-Atherton model of four hysteresis loops with increase of magnetizing field amplitude
% 
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Biedrzycki R., Jackiewicz D., Szewczyk R. "Reliability and Efficiency of Differential Evolution Based Method of Determination 
%     of Jiles-Atherton Model Parameters for X30Cr13 Corrosion Resisting Martensitic Steel" Journal of Automation, Mobile Robotics and Intelligent Systems 8 (2014) 63-68.
%
% USAGE:
% demo05_octave_view
% 
% IMPORTANT: Demo requires "odepkg" package installed and loaded  
%


clear all
clc

page_screen_output(0);
page_output_immediately(1);  % print immediately at the screen

fprintf('\n\nView of identified parameters of Jiles-Atherton models parameters for four hysteresis loops.');
fprintf('\nDemonstration optimized for OCTAVE. For MATLAB please use demo05_matlab_view.m ');
fprintf('\nDemonstration requires odepkg, struct and optim packages installed.\n\n');

% check if odepkg is installed. Load odepkg if installed, but not loaded.
ChkPkg('odepkg');

% check if struct is installed. Load odepkg if installed, but not loaded.
ChkPkg('struct');

% check if optim is installed. Load odepkg if installed, but not loaded.
ChkPkg('optim');

% Load measured B(H) characterisitcs of Mn-Zn ferrite

cd ('Characterisitcs_isotropic_mat');
load('H_MnZn_ferrite.mat');
load('B_MnZn_ferrite.mat');
cd ('..');
 
fprintf('Load measured B(H) characterisitcs of Mn-Zn ferrite... done\n\n');

% prepare data 

mi0=4.*pi.*1e-7;

load('demo05_results.mat');

func = @(JApointn) JAn_loops_target( JApointn, JApoint0, HmeasT, BmeasT, 1);

fprintf('\n\nCalculations...\n\n');

Ftarget=func(JApoint_res);

BsimT = JAn_loops(JApoint0(1).*JApoint_res(1), JApoint0(2).*JApoint_res(2), JApoint0(3).*JApoint_res(3), ...
 JApoint0(4).*JApoint_res(4), JApoint0(5).*JApoint_res(5), HmeasT, 1 );

fprintf('Results of optimisation:\n'); 
fprintf('Target function value: Ftarget=%f\n',Ftarget);
fprintf('JA model params: a=%f(A/m), k=%f(A/m), c=%f, Ms=%e(A/m), alpha=%e \n\n',  ...
 JApoint0(1).*JApoint_res(1), JApoint0(2).*JApoint_res(2), JApoint0(3).*JApoint_res(3), ...
 JApoint0(4).*JApoint_res(4), JApoint0(5).*JApoint_res(5));
 
plot(HmeasT, BmeasT,'r',HmeasT,BsimT,'k');
xlabel('H (A/m)');
ylabel('B (T)');
grid;

