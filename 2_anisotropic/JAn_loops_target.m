function Ftarget = JAn_loops_target(JApointn,JApoint0,HmeasT,BmeasT,SolverType,FixedStep)
%
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
% JAn_loop_target_mod calculates target function for optimisation. 
% Target function is the sum of squared differences between results of modelling and results of measurements.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C., Atherton D. "Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Szewczyk R. "Computational problems connected with Jiles-Atherton model of magnetic hysteresis". Advances in Intelligent Systems and Computing (Springer) 267 (2014) 275.
% USAGE:
% [H,M] = JAsolver(a,k,c,Ms,alpha,Hstart,Hend,M0,SolverType)
% 
% IMPORTANT: H must be a COLUMN vector 
%
% INPUT:
% JApoint - point for optimisation [a k c Ms alpha Kan psi] (vector)
% HmeasT  - results of measurements of magnetizing field, A/m (matrix - set of data in columns)
% Bmeast  - results of measurements of flux density, T (matrix - set of data in columns)
% SolverType - select the solver for ODE
%               1 - ode23()
%               2 - ode45()
%               3 - ode23s()
%               4 - rk4() 
% FixedStep - select the solver for integration:
%               1: quadtrapz(), 0: quadgk()
%
% OUTPUT:
% Ftarget -  sum of squared differences between results of measurements BmeasT and 
%            results of simulation BsimT - target function for optimisation (scalar)


if size(HmeasT)~=size(BmeasT)
    fprintf('\n\n*** ERROR in JAn_loops_target: diferent sizes HmeasT and BmeasT.\n\n');
    Ftarget=1e30;
    return
end  

JApoint=JApointn.*JApoint0;

if numel(JApoint)~=7
    fprintf('\n\n*** ERROR in JAn_loops_target: wrong size of JApoint.\n\n');
    Ftarget=1e30;
    return
end  

a=JApoint(1);
k=JApoint(2);
c=JApoint(3);
Ms=JApoint(4);
alpha=JApoint(5);
Kan=JApoint(6);
psi=JApoint(7);

if ((a<0) || (k<0) || (c<0) || (c>1) || (alpha<0) || (Kan<0)) 
    fprintf('\n\n*** ERROR in JAn_loops_target: unphysical parameters.\n\n');
    Ftarget=1e30;
    return
end           % check if parameters are physical

BsimT = JAn_loops(a,k,c,Ms,alpha,Kan,psi,HmeasT, SolverType, FixedStep);

Ftarget=sum(sum((BmeasT-BsimT).^2));  % calculate target function

if isnan(Ftarget)==1
    Ftarget=1e30;                   % check if the target function is NaN. In such a case make it very high
end 

end
