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
% 
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% 
% USAGE:
% demo01_matlab_simple_loop
% 


clear all
clc

fprintf('\n\nSimple demonstration of Jiles-Atherton model.');
fprintf('\nDemonstration optimized for MATLAB. For OCTAVE please use demo01_octave_simple_loop.m\n\n');


% Calculate Jiles-Atheton model for the case presented in [1]
 
% prepare variables for modelling

Ms=1.6e6;
a=1100;
alpha=1.6e-3;
k=400;
c=0.2;        % Parameters of Jiles-Atherton model specified in [1]

fprintf('Calculation for parameters: \na=%f(A/m), k=%f(A/m), c=%f, Ms=%e(A/m), alpha=%e \n\n',a,k,c,Ms,alpha);

H=[0:60:6000 6000:-60:-6000 -6000:60:6000]'; % magnetizing field H - column vector

M0=0;         % sample demagnetized at the beginning
SolverType=1; % ode23() solver

fprintf('Calculations started. Expected time of calculations: below 1 minute.\n');
[Hmodel,Bmodel] = JAsingle_loop(a,k,c,Ms,alpha,H,M0,SolverType);
fprintf('Calculations finished.\n\n');

plot(Hmodel, Bmodel,'-o');
xlabel('H (A/m)');
ylabel('B (T)');
grid;


