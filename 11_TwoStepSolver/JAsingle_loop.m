function [Hw,Bw] = JAsingle_loop(Hmeas,a,k,c,Ms,alpha,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType)
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
% Hmeas - magnetizing field, A/m (the column vector, A/m)
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
%
% OUTPUT:
% Hw - set of output values of magnetizing field H, A/m (vector)
% Bw - set of output values of flux density B, T (vector)

mi0=4.*pi.*1e-7;
M0=0;
H=Hmeas;


% determine default values for not determined variables
if nargin<13 || isempty(IntType)
   IntType=0;
end
if nargin<12 || isempty(psi)
   psi=0;
end
if nargin<11 || isempty(Kan)
   Kan=0;
end
if nargin<10 || isempty(AnisoType)
   AnisoType=0;
end
if nargin<9 || isempty(IsoAniso)
   IsoAniso=0;
end
if nargin<8 || isempty(SolverType)
   SolverType=0;
end
if nargin<7 || isempty(ModelType)
   ModelType=0;
end


if ((a<0) || (k<0) || (c<0) || (c>1) || (alpha<0) || (Kan<0)) 
    fprintf('\n\n*** ERROR in JAsingle_loops: unphysical parameters.\n\n');
    Hw=H;
    Bw=zeros(size(H));
    return
end           % check if parameters are physical


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
   [Hi, Mi]=JAsolver(a,k,c,Ms,alpha,H(ip),H(ik),M0,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType);
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
         [Hi, Mi]=JAsolver(a,k,c,Ms,alpha,H(ip),H(ik),M0,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType);
         Hi(isnan(Hi))=1;
         Mi(isnan(Mi))=1;
         if numel(Hi)>2
            Mw(ip+1:ik)=interp1(Hi,Mi,H(ip+1:ik),'cubic','extrap');
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
         [Hi, Mi]=JAsolver(a,k,c,Ms,alpha,H(ip),H(ik),M0,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType);
         Hi(isnan(Hi))=1;
         Mi(isnan(Mi))=1;
         if numel(Hi)>2
            Mw(ip+1:ik)=interp1(Hi,Mi,H(ip+1:ik),'cubic','extrap');
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

[Hi, Mi]=JAsolver(a,k,c,Ms,alpha,H(ip),H(ik),M0,ModelType,SolverType,IsoAniso,AnisoType,Kan,psi,IntType);
Hi(isnan(Hi))=1;
Mi(isnan(Mi))=1;
if numel(Hi)>2
   Mw(ip+1:ik)=interp1(Hi,Mi,H(ip+1:ik),'cubic','extrap');
else
   Mw(ip+1)=Mi(end);
end

Hw(ip+1:ik)=H(ip+1:ik);
M0=Mi(end);


Bw=Mw.*mi0+Hw.*mi0;

end

