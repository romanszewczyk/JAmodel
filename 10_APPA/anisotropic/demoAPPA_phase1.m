% The MIT License (MIT)
%
% Copyright (c) 2016 Roman Szewczyk, Peng Cheng
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
% Demonstration of identification of parameters of different types of Jiles-Atherton model of four hysteresis loops with increase of magnetizing field amplitude
% fminsearch() function used with Nelder and Mead Simplex algorithm (a derivative-free method)
%
% HmeasT - magnetizing field, A/m (set of column vectors - matrix, A/m)
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
% AnisoType - select the type of anisotropy model
%               0: uniaxial anisotropy
%               1: GO anisotropy
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
% AUTHORS: Roman Szewczyk, rszewczyk@onet.pl, Peng Cheng, Peng.Cheng@miun.se
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Szewczyk R. "Computational problems connected with Jiles-Atherton model of magnetic hysteresis". Advances in Intelligent Systems and Computing (Springer) 267 (2014) 275.
%
% IMPORTANT: Demo requires "odepkg", "struct" and "optim" packages installed and loaded  
%

clear all
clc

page_screen_output(0);
page_output_immediately(1);  % print immediately at the screen


fprintf('\n\nDemonstration of identification of different types of Jiles-Atherton models parameters for four hysteresis loops.');
fprintf('\nDemonstration optimized for OCTAVE. ');
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

Ms0=1.5.*max(max(BmeasT))./mi0;
a0=20;
alpha0=1e-6;
k0=12;
c0=0.8; 
Kan=250;
psi=0;
       % Initial parameters of Jiles-Atherton model for optimisation

JApoint0=[a0 k0 c0 Ms0 alpha0 Kan];

ModelType=1;
SolverType=1;
IsoAniso=1;
AnisoType=1;
IntType=0;

func = @(JApointn) JAn_loops_target( JApointn, JApoint0, HmeasT, BmeasT, ModelType,SolverType,IsoAniso,AnisoType,IntType);
options=optimset('Display','iter','MaxFunEvals',500);

fprintf('Optimization process started... (first cycle expected in less than 10 min.)\n\n');

tic

JApoint_res=fminsearch(func,[1 1 1 1 1 1],options);

toc

fprintf('\n\nOptimiation process done.\n\n');

Ftarget=func(JApoint_res);

BsimT = JAn_loops(HmeasT,JApoint0(1).*JApoint_res(1), JApoint0(2).*JApoint_res(2), JApoint0(3).*JApoint_res(3), ...
 JApoint0(4).*JApoint_res(4), JApoint0(5).*JApoint_res(5), ModelType,SolverType,IsoAniso,AnisoType, JApoint0(6).*JApoint_res(6),psi,IntType);
 

fprintf('Results of optimisation:\n'); 
fprintf('Target function value: Ftarget=%f\n',Ftarget);
fprintf('JA model params: a=%f(A/m), k=%f(A/m), c=%f, Ms=%e(A/m), alpha=%e, Kan=%e(J/m3), psi=0 \n\n',  ...
 JApoint0(1).*JApoint_res(1), JApoint0(2).*JApoint_res(2), JApoint0(3).*JApoint_res(3), ...
 JApoint0(4).*JApoint_res(4), JApoint0(5).*JApoint_res(5), JApoint0(6).*JApoint_res(6));
 
fprintf('Optimisation done.\n\n');

plot(HmeasT, BmeasT,'or',HmeasT,BsimT,'k');
xlabel('H (A/m)');
ylabel('B (T)');
grid;

JApoint_optim=JApoint0.*JApoint_res;

save -v7 demo03_results.mat JApoint_optim JApoint0 JApoint_res

