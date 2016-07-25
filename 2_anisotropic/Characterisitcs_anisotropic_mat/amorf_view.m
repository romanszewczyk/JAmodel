clear all

load('H_amorf.mat');
load('B_amorf.mat');

plot(HmeasT,BmeasT);
title('Magnetic hysteresis loops B(H) of M-680 amorphous alloy core with strong perpendicular anisotropy');
xlabel('H (A/m)');
ylabel('B (T)');
grid;
