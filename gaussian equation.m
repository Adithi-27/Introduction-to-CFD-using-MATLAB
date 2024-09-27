clc
clear all
format compact
close all

%variables
x = 0:1:20;
mu = 10;
sigma = 3;

%gaussian equation
y = 1/(sigma*sqrt(2*pi))*exp(-((x-mu).^2)/(2*sigma^2));

%outputs
y(1);
plot(x,y)