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
% Demonstration of identification of parameters of Jiles-Atherton model of four hysteresis loops with increase of magnetizing field amplitude
% 
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Szewczyk R. "Computational problems connected with Jiles-Atherton model of magnetic hysteresis". Advances in Intelligent Systems and Computing (Springer) 267 (2014) 275.
%
% USAGE:
% demo04_octave_view
% 
% IMPORTANT: Demo requires "odepkg" package installed and loaded  
%


clear all
clc

page_screen_output(0);
page_output_immediately(1);  % print immediately at the screen


fprintf('\n\nDemonstration of identification of Jiles-Atherton models parameters for four hysteresis loops.');
fprintf('\nDemonstration optimized for OCTAVE. For MATLAB please use demo04_octave_view.m ');
fprintf('\nDemonstration requires odepkg, struct and optim packages installed.\n\n');


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

% check if struct is installed. Load odepkg if installed, but not loaded.

inst_pkg = pkg ("list");

[i,j]=size(inst_pkg);

struct_inst=0;
struct_loaded=0;

for i=1:j
    if size(findstr(inst_pkg{1,i}.name,'struct'))>0
       struct_inst=1;
       if inst_pkg{1,i}.loaded==1;
          struct_loaded=1;
       end
    end
end

if struct_inst==0
   fprintf('\n *** ERROR: struct must be installed to solve ODEs.\n To solve problem try: pkg install -forge struct\n\n');
   return
else
   fprintf('\n struct installed...ok.');
end
   
 if struct_loaded==0
   fprintf('\n WARNING: struct is installed but not loaded.\n');
   pkg load struct
   fprintf(' Problem solved: struct is loaded now.\n\n');
   else
   fprintf('\n struct loaded...ok.\n\n');
end

% check if optim is installed. Load odepkg if installed, but not loaded.

inst_pkg = pkg ("list");

[i,j]=size(inst_pkg);

optim_inst=0;
optim_loaded=0;

for i=1:j
    if size(findstr(inst_pkg{1,i}.name,'optim'))>0
       optim_inst=1;
       if inst_pkg{1,i}.loaded==1;
          optim_loaded=1;
       end
    end
end

if optim_inst==0
   fprintf('\n *** ERROR: optim must be installed to solve ODEs.\n To solve problem try: pkg install -forge odepkg\n\n');
   return
else
   fprintf('\n optim installed...ok.');
end
   
 if optim_loaded==0
   fprintf('\n WARNING: optim is installed but not loaded.\n');
   pkg load optim
   fprintf(' Problem solved: optim is loaded now.\n\n');
   else
   fprintf('\n optim loaded...ok.\n\n');
end

% Load measured B(H) characterisitcs of Mn-Zn ferrite

cd ('Characterisitcs_isotropic_mat');
load('H_MnZn_ferrite.mat');
load('B_MnZn_ferrite.mat');
cd ('..');
 
fprintf('Load measured B(H) characterisitcs of Mn-Zn ferrite... done\n\n');

% prepare data 

mi0=4.*pi.*1e-7;

load('demo04_results.mat');

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

