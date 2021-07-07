%% Nonlinear Cobweb
% This code calculates the following for a simple non-linear cobweb model: 
% 1. Iteration Loop
% 2. Time Series 
% 3. Bifurcation Diagram
% 4. Parameter Basins of Attraction 
% This is inspired by Hommes (2013) and his lectures at university of
% Amsterdam as well as the UvA CeNDEF software 'E&FChaos'. 
% (Written by Kian Abbas Nejad)

%% The Model 
% This 1D Model is given by the following supply and demand schedules: 
% $D(p_t)=a-b p_t$
% $S(p_t^e)=\arctan [\lambda(p_t^e-c)]+d$
% And expectations are specified adaptively by: 
% p_t^e = (1-w)p^e_{t-1} + w p_{t-1}
% Which result in the following law of motion: 
%
% $p^e_{t+1}=f(p_t^e)=(1-w)p_t^e+w \frac{a-d-\arctan[\lambda(p_t^e-c)]}{b}$
%% Housekeeping

clearvars
close all
clc

% Parameters
w = 0.5;
b = 0.25;
lambda = 4.4;
a = 1;
c = 0;
d = 0;

% Settings
Pmin = -10;
Pmax = 10;
step = 0.01;
n = 3; %number of iterations
p0 = -6; %initial condition
p1 = p0; %vector for iteration loop
num = 5; % number of iterations for plotting
P = Pmin:step:Pmax;
I = P; % vector for iterations

% The Map
fp = @(p) (1-w).*p + w.*(a-d-atan(lambda.*(p-c)))./b;

%% Time Series

tmax = 50;
trans = 0; %number of transient iterations not to plot
T = 1:tmax;
P = zeros(tmax,1);

P(1)= p0; % init. cond. is p0 variable

for i=2:tmax
    pp = P(i-1);
    P(i)=fp(pp);
end

fig2 = figure(2);
plot(T(trans+1:tmax),P(trans+1:tmax),'k')

fig2.Name= 'Time Series';
ylabel('$p_t$','FontSize',16,'interpreter','latex');
xlabel('$t$','FontSize',16,'interpreter','latex');
title('Time Series of $p_t$','FontSize',14,'interpreter','latex');

%% Iteration Loop

IterationMapFun(fp,n,p0,num,'p',[-10 10 -10 10]) 

%% Bifurcation Diagram

fp = @(p,W) (1-W).*p + W.*(a-d-atan(lambda.*(p-c)))./b;

BifurcationFun(fp,0.2:0.001:0.9,500,1000,'p','w')

%% Basin of Attraction

Fp = @(p,L,W) (1-W).*p + W.*(a-d-atan(L.*(p-c)))./b;

BasinFun(Fp,7,0:0.002:10,0.1:0.002:0.8,6,300,400,'\lambda','w')
