%HW3, Q2
clear all; close all; clc

x_begin= 1.920;
x_end=2.080;
x_step=0.001;

x=x_begin:x_step:x_end;
y1=zeros(size(x));
y2=zeros(size(x));


p=[1 -18 144 -672 2016 -4032 5376 -4608 2304 -512] % polynomial function
y1 = polyval(p,x)
y2 = (x-2).^9;

plot(x,polyval(p,x),'b')
hold on
plot(x,y2,'r')
xlabel('x','FontSize', 20)
ylabel('p(x)','FontSize', 20)
title('Polynominal p(x)','FontSize', 20)
legend({'right-hand side','left-hand side'},'Location','northwest')