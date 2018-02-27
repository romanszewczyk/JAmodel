% The MIT License (MIT)
%
% Copyright (c) 2017 Roman Szewczyk
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
% Demonstration of two stage identification of parameters of Jiles-Atherton model 
% with the hysteresis described like in original model and with Cheng et al. modification
% 
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] P. Cheng, R. Szewczyk "Modified description of magnetic hysteresis in Jiles-Atherton model"
%
% Paper presented during the Automation2018 Conference, automation.piap.pl
% 
%
% IMPORTANT: Demo requires "odepkg", "struct" and "optim" packages installed and loaded  
%


clear all
diary on
clc

page_screen_output(0);
page_output_immediately(1);  % print immediately at the screen


fprintf('\n\nDemonstration of identification of Jiles-Atherton models parameters for three hysteresis loops.');
fprintf('\nMultistep solver.');
fprintf('\nDemonstration optimized for OCTAVE.  ');
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

% Initial analyse of input data
fprintf('Initial analyse of input data... \n\n');


mi0=4.*pi.*1e-7;


Number_of_loops=size(HmeasT,2);       % number of loops

[Hmax_colV Hmax_colI] = max(HmeasT);
[Hmax_V Hmax_rowI] = max(Hmax_colV);
Hmax_colI=Hmax_colI(Hmax_rowI);


Hmax_V=HmeasT(Hmax_colI,Hmax_rowI);   % test - find the maximal H value
Bmax_V=BmeasT(Hmax_colI,Hmax_rowI);   % find the maximal B value
Mmax_V=Bmax_V./mi0;                   % calculate maximal M value


HcT=[];
Hmax_loop=HmeasT(:,Hmax_rowI);
Bmax_loop=BmeasT(:,Hmax_rowI);


for i=2:numel(Hmax_loop)
    
    x1=Hmax_loop(i-1);
    x2=Hmax_loop(i);
    
    y1=Bmax_loop(i-1);
    y2=Bmax_loop(i);
    
    if (((y1>0) && (y2<0)) || ((y1<0) && (y2>0)))  
       
       a=(y2-y1)./(x2-x1);
       b=y1-a.*x1;
       Hc_i=-1.*b./a;
       
       HcT=[HcT Hc_i];
       
       end
   end
   
Hc= mean(abs(HcT));

 
fprintf('Number of loops: %i \n',Number_of_loops);
fprintf('Maximal magnetizing field H(A/m): %5.2f, pos. at col:%i , pos. at row:%i \n',Hmax_V, Hmax_colI,Hmax_rowI);
fprintf('Maximal magnetization M(A/m): %6.0f \n',Mmax_V);
fprintf('Maximal flux density B(T): %1.3f \n',Bmax_V); 
fprintf('Number of x axis crossings at major loop: %i \n',numel(HcT)); 
fprintf('Estimated value of Hc (A/m): %3.2f \n\nDone. \n\n',Hc); 

figure(1)
plot(HmeasT,BmeasT);
xlabel('H (A/m)');
ylabel('B (T)');
title('Input data');
grid;
drawnow;


% Identification of parameters for anhysteretic loop
fprintf('Identification of parameters for anhysteretic loop...\n\n');

a=Hc;
Ms=Mmax_V;
alpha=1e-6;

MahT = Mah_loop(Hmax_loop,a,Ms,alpha); 

BahT=(MahT+Hmax_loop).*mi0;

figure(2)
plot(Hmax_loop,Bmax_loop,Hmax_loop,BahT);
xlabel('H (A/m)');
ylabel('B (T)');
title('Anhysteretic curve starting point');
grid;
drawnow;


func=@(opt) Mah_loop_target(Hmax_loop,Bmax_loop,opt(1),opt(2),opt(3));

fprintf('Initial target value for a=%2.3f (A/m), Ms=%6.0f (A/m), alpha=%e Target=%2.3f \n\n',a,Ms,alpha, func([a Ms alpha])); 

ctl.XVmin = [0.2.*a Ms 1e-8];
ctl.XVmax = [5.*a 2.*Ms 5e-4];
ctl.refresh = 1;
ctl.maxiter = 25;
ctl.constr = 1;
ctl.NP = 100;


fprintf('Optimization process started... \n\n');

tic

[JAMah_res, obj_value, nfeval, convergence] = de_min (func, ctl);

toc

fprintf('\n\nOptimiation process done.\n\n');

fprintf('Final target value for a=%2.3f (A/m), Ms=%6.0f (A/m), alpha=%e Target=%2.3f \n\n',JAMah_res(1), JAMah_res(2), JAMah_res(3), func(JAMah_res)); 

MahT2 = Mah_loop(Hmax_loop,JAMah_res(1), JAMah_res(2), JAMah_res(3));

BahT2=(MahT2+Hmax_loop).*mi0;

figure(3)
plot(Hmax_loop,Bmax_loop,Hmax_loop,BahT2);
xlabel('H (A/m)');
ylabel('B (T)');
title('Anhysteretic curve optimised');
grid;
drawnow;

% Adding the hysteresis

fprintf('Adding the hysteresis: JA model... \n\n');
fprintf('Initial point... \n');

k=Hc;
c=0.7;

ModelType=0;
SolverType=4;


BsimT3 = JAn_loops(HmeasT,JAMah_res(1),k,c,JAMah_res(2),JAMah_res(3));

figure(4)
plot(HmeasT,BmeasT,HmeasT,BsimT3);
xlabel('H (A/m)');
ylabel('B (T)');
title('Hysteretic curve initial');
grid;
drawnow;


func=@(kc) JAn_loops_target([1 1 1 1 1 1 1],[JAMah_res(1) kc(1) kc(2) JAMah_res(2) JAMah_res(3) 0 0],HmeasT,BmeasT,ModelType,SolverType);

ctl.XVmin = [0.2.*a 1e-6];
ctl.XVmax = [5.*a 1];
ctl.refresh = 1;
ctl.maxiter = 25;
ctl.constr = 1;
ctl.NP = 30;

fprintf('Optimization process started... \n\n');

tic

[JA_res, obj_value, nfeval, convergence] = de_min (func, ctl);

toc

fprintf('\n\nOptimiation process done.\n\n');

BsimT4 = JAn_loops(HmeasT,JAMah_res(1),JA_res(1),JA_res(2),JAMah_res(2),JAMah_res(3),ModelType,SolverType);

figure(5)
plot(HmeasT,BmeasT,HmeasT,BsimT4);
xlabel('H (A/m)');
ylabel('B (T)');
title('Hysteretic curve after optimisation of k and c');
grid;
drawnow;


clear func 
save -v7 ResultsJA.mat *

% Adding the hysteresis

fprintf('Adding the hysteresis: Cheng model... \n\n');
fprintf('Initial point... \n');

k=Hc;
c=0.7;

ModelType=3;
SolverType=4;


BsimT3 = JAn_loops(HmeasT,JAMah_res(1),k,c,JAMah_res(2),JAMah_res(3));

figure(6)
plot(HmeasT,BmeasT,HmeasT,BsimT3);
xlabel('H (A/m)');
ylabel('B (T)');
title('Hysteretic curve initial');
grid;
drawnow;


func=@(kc) JAn_loops_target([1 1 1 1 1 1 1],[JAMah_res(1) kc(1) kc(2) JAMah_res(2) JAMah_res(3) 0 0],HmeasT,BmeasT,ModelType,SolverType);

ctl.XVmin = [0.2.*a 1e-6];
ctl.XVmax = [5.*a 1];
ctl.refresh = 1;
ctl.maxiter = 25;
ctl.constr = 1;
ctl.NP = 30;

fprintf('Optimization process started... \n\n');

tic

[JA_res, obj_value, nfeval, convergence] = de_min (func, ctl);

toc

fprintf('\n\nOptimiation process done.\n\n');

BsimT4 = JAn_loops(HmeasT,JAMah_res(1),JA_res(1),JA_res(2),JAMah_res(2),JAMah_res(3),ModelType,SolverType);

figure(7)
plot(HmeasT,BmeasT,HmeasT,BsimT4);
xlabel('H (A/m)');
ylabel('B (T)');
title('Hysteretic curve after optimisation of k and c');
grid;
drawnow;


clear func 
save -v7 ResultsCheng.mat *

diary off
