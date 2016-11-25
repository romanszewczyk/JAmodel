function [Hw,Bw] = JAsingle_loop(a,k,c,Ms,alpha,Kan,psi,t,H,M0_,SolverType,FixedStep)
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
% JAsingle_loop splits the given magnetizing field H into monotonous sections and control solving
% ODE stating Jiles-Atherton model
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
% a    - quantifies domain density, A/m (scalar)
% k    - quantifies average energy required to break pinning side, A/m (scalar)
% c    - magnetization reversibility, 0..1 (scalar)
% Ms   - saturation magnetization, A/m (scalar)
% alpha - Bloch coefficient (scalar)
% H    - magnetizing field, A/m (column vector)
% M0   - initial value of magnetization for ODE, A/m (scalar)
% SolverType - select the solver for ODE
%               1 - ode23()
%               2 - ode45()
%               3 - ode23s()
%               4 - rk4() 
%
% OUTPUT:
% Hw - set of output values of magnetizing field H, A/m (vector)
% Bw - set of output values of flux density B, T (vector)

mi0=4.*pi.*1e-7;
M0=M0_;

if size(H,2)>1
   fprintf('\n\n***ERROR in JAsingle_loop: H must be the column vector\n\n');
   Bw=0;
   Hw=0;
   return
end

if numel(H)<2
   fprintf('\n\n***ERROR in JAsingle_loop: H must be at least two elements vector\n\n');
   Bw=0;
   Hw=0;
   return
end

if numel(H)==2
   [Hi, Mi]=JAsolver(a,k,c,Ms,alpha,Kan,psi,t,H(ip),H(ik),M0,SolverType,FixedStep);
   Hi(isnan(Hi))=1;
   Mi(isnan(Mi))=1;
   Hw=H;
   Mw=[M0 Mi(end)];
   Bw=Mw.*mi0+Hw.*mi0;
   return
   end


ip=1;
ik=2;

Hw=zeros(size(H));
Hw(1)=H(1);

Mw=zeros(size(H));
Mw(1)=M0;


while ik<numel(H)

   if H(ik)>H(ip)
      if H(ik+1)>=H(ik)
         ik=ik+1;
      else
         [Hi, Mi]=JAsolver(a,k,c,Ms,alpha,Kan,psi,t,H(ip),H(ik),M0,SolverType,FixedStep);
         Hi(isnan(Hi))=1;
         Mi(isnan(Mi))=1;
         if numel(Hi)>2
            Mw(ip+1:ik)=interp1(Hi,Mi,H(ip+1:ik),'cubic');
         else
            Mw(ip+1)=Mi(end);
         end
         Hw(ip+1:ik)=H(ip+1:ik);
         M0=Mi(end);

         ip=ik;
         ik=ip+1;
      end
   end

   if H(ik)<H(ip)
      if H(ik+1)<=H(ik)
         ik=ik+1;
      else
         [Hi, Mi]=JAsolver(a,k,c,Ms,alpha,Kan,psi,t,H(ip),H(ik),M0,SolverType,FixedStep);
         Hi(isnan(Hi))=1;
         Mi(isnan(Mi))=1;
         if numel(Hi)>2
            Mw(ip+1:ik)=interp1(Hi,Mi,H(ip+1:ik),'cubic');
         else
            Mw(ip+1)=Mi(end);
         end
         Hw(ip+1:ik)=H(ip+1:ik);
         M0=Mi(end);
         
         ip=ik;
         ik=ip+1;
      end
   end
   
   if H(ik)==H(ip)
      ik=ik+1;
   end

end

[Hi, Mi]=JAsolver(a,k,c,Ms,alpha,Kan,psi,t,H(ip),H(ik),M0,SolverType,FixedStep);
Hi(isnan(Hi))=1;
Mi(isnan(Mi))=1;
if numel(Hi)>2
   Mw(ip+1:ik)=interp1(Hi,Mi,H(ip+1:ik),'cubic');
else
   Mw(ip+1)=Mi(end);
end

Hw(ip+1:ik)=H(ip+1:ik);
M0=Mi(end);


Bw=Mw.*mi0+Hw.*mi0;

end

