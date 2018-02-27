function [t_res, y_res]=rk4(fun,time,f0,N)
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
% rk4 is the function for fixed steps solving ODEs by Runge-Kutta 4-th order method
% WARNING: solver is fast, but accuracy is not controlled
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
%
% USAGE:
% [t_res, y_res]=rk4(fun,time,f0,N)
% 
% INPUT:
% fun    - function for ODE
% time   - time for solver: [time_start time_stop] (vector)
% f0     - value of the solution for time_start, (scalar)
% N      - fixed number of steps, (scalar)
%
% OUTPUT:
% t_res  - time points for solutions (termined by time_start, time_stop and N) (vector)
% y_res  - solution of the ODE for in t_res points
%

step=(time(2)-time(1))./N;
t_res=time(1):step:time(2);
y_res=t_res.*0;

y_res(1)=f0;

for i=1:N
   a1=step.*fun(t_res(i),y_res(i));
   a2=step.*fun(t_res(i)+step./2,y_res(i)+a1./2);
   a3=step.*fun(t_res(i)+step./2,y_res(i)+a2./2);
   a4=step.*fun(t_res(i)+step,y_res(i)+a3);
   y_res(i+1)=y_res(i)+(a1+2.*a2+2.*a3+a4)./6;
end


end  % end of function 
