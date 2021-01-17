function res=mu_i_GO(Ms, a, Kan, psi)
%
% The MIT License (MIT)
%
% Copyright (c) 2020 Roman Szewczyk
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
% function for calculation the anhysteretic curve of materials with GO anisotropy
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% [1] Jiles D. C. , Atherton D. "Theory of ferromagnetic hysteresis" Journal of Magnetism and Magnetic Materials 61 (1986) 48.
% [2] Ramesh A., Jiles D., Roderik J. "A model of anisotropic anhysteretic magnetization" IEEE Transactions on Magnetics 32 (1996) 4234–4236. [Google Scholar] [CrossRef]
% [3] Ramesh A., Jiles D., Bi Y. "Generalization of hysteresis modeling to anisotropic materials" Journal of Applied Physics  81 (1997) 5585–5587.
% [4] Szewczyk R. "Validation of the Anhysteretic Magnetization Model for Soft Magnetic Materials with Perpendicular Anisotropy" Materials  7 (2014) 5109-5116.
%


mi0=4.*pi.*1e-7;    % magnetic constant

H=0.001;            % (A/m) - numerical limit to 0 A/m


if ((Ms>0) && (a>0) && (Kan>0))
 
    E1= @(t) H./a.*cos(t)-Kan/(Ms.*mi0.*a).*((cos(psi-t)).^2.*(sin(psi-t)).^2+0.25.*(sin(psi-t)).^4);
    E2= @(t) H./a.*cos(t)-Kan/(Ms.*mi0.*a).*((cos(psi+t)).^2.*(sin(psi+t)).^2+0.25.*(sin(psi+t)).^4);

    F1= @(t) exp((E1(t)+E2(t))./2).*sin(t).*cos(t);
    F2= @(t) exp((E1(t)+E2(t))./2).*sin(t);

    res = Ms./H.*quad(F1,0,pi)./quad(F2,0,pi);
else
    res=0;
end
    

end
