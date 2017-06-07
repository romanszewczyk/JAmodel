clear all

load('H_M130_27s.mat');
load('B_M130_27s.mat');

plot(HmeasT,BmeasT);
title('Magnetic hysteresis loops B(H) of M130-27s electric steel measured in the easy axis direction');
xlabel('H (A/m)');
ylabel('B (T)');
grid;
