clear all

load('H_MnZn_ferrite.mat');
load('B_MnZn_ferrite.mat');

plot(HmeasT,BmeasT);
title('Magnetic hysteresis loops B(H) of Mn-Zn ferrite');
xlabel('H (A/m)');
ylabel('B (T)');
grid;
