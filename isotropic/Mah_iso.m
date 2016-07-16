function result = Mah_iso(a,Ms,He)
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
% Mah_iso is the function for calculation of anhysteretic, isotropic magnetization in Jiles-Atherton model
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% Jiles D. C. , Atherton D. „Theory of ferromagnetic hysteresis” Journal of Applied Physics 55 (1984) 2115.
% Jiles D. C. , Atherton D. „Theory of ferromagnetic hysteresis” Journal of Magnetism and Magnetic Materials 61 (1986) 48.
%
% USAGE:
% M = Mah_iso(a, Ms, He)
% 
% INPUT:
% a    - quantifies domain density, A/m (scalar)
% Ms   - saturation magnetization, A/m (scalar)
% He   - effective magnetizing field, A/m (scalar, vector or matrix)
% 
% OUTPUT:
% result - value of anhysteretic, isotropic magnetization in Jiles-Atherton model, A/m (scalar, vector or matrix, size of He)
%

if (max(size(a))>1 || max(size(Ms)>1))
  
  fprintf('\n*** ERROR in Mah_iso: a or Ms is not a scalar. ***\n');
  return;
  
end
  
if (Ms==0) || (a==0)

  result=zeros(size(He));
  return;

end  % of if


i=(He==0);        % find singularities for He==0

He(i)=1;          % remove singularities

result=Ms.*(coth(He./a)-(a./He));   % calculate anhysteretic, isotropic magnetization

result(i)=0;      % correct values for singularities


end  % of function

