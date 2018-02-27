function [xopt, func_val] = GradDescSolver(f,x0,theta,MaxIter,MinGrad)
%
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
% Gradient descent solver
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% https://en.wikipedia.org/wiki/Trapezoidal_rule
%
% USAGE:
% res = quadtrapz(f,xstart,xend,N)
% 
% INPUT:
% f       - optimised function
% x0  - starting point of integration 
%
% OUTPUT:
% xopt
%

fprintf('\n\nGradient descent optimisation solver');
fprintf('\nStart value: %e\n\n',f(x0));

fold=1e99;
fnew=f(x0);

g=ones(size(x0));
x=x0;

iter=0;

while ((iter<MaxIter) && (fnew<fold) && (max(max(g))>MinGrad))

% calculate gradients

for i=1:numel(x0)
     
     x1=x;
     x2=x;
     
     x1(i)=x1(i)+1e-8.*x1(i);
     x2(i)=x2(i)-1e-8.*x2(i);
     
     g(i)=(f(x1)-f(x2))./(2e-8.*x1(i));
     
     if g(i)>100
        g(i)=100;
        end
     if g(i)<-100
        g(i)=-100;
        end
        
        
end

xopt=x+x.*theta.*g;

fold=fnew;
fnew=f(xopt);

if fnew<fold
  x=xopt;
  else
  xopt=x;
  end

iter=iter+1;

fprintf('\nIter: %i, value: %e, gmax: %e \n',iter ,fnew ,max(max(g)));

% display(x);

end

func_val=fnew;

endfunction
