%% Nonlinear Cobweb
% This code calculates the following for a simple non-linear cobweb model: 
% 1. Iteration Loop
% 2. Time Series 
% 3. Bifurcation Diagram
% 4. Parameter Basins of Attraction (incomplete)
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
p0 = 4.75; %initial condition
p1 = p0; %vector for iteration loop
num = 9; % number of iterations for plotting
P = Pmin:step:Pmax;
I = P; % vector for iterations

% The Map
fp = @(p) (1-w).*p + w.*(a-d-atan(lambda.*(p-c)))./b;

%% Iteration Loop
fig1= figure(1);
hold on
grid on
fig1.Name = 'Iteration';
ylabel('$p_{t+1}$','FontSize',16,'interpreter','latex')
xlabel('$p_t$','FontSize',16,'interpreter','latex')
title('Iteration of $p_{t+3}=f^3(p_t)$','FontSize',14,'interpreter','latex')
axis([Pmin Pmax Pmin Pmax])
plot(P,P,'b--')

for i=1:n
    I=fp(I);
end

% Plotting the nth iteration of the map
plot(P,I,'r')

for i=1:num
    for k=1:n
        p1 =fp(p1);
    end
    plot([p0 p0],[p0 p1],'k')
    plot([p0 p1],[p1 p1],'k')
    p0 = p1;
end

%% Time Series
tmax = 100;
trans = 50; %number of transient iterations not to plot
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

%% Bifurcation Diagram

pmin = 0.2; %parameter minimum
pmax = 0.9; %parameter max
pstep = 0.001; %parameter step size
tmax = 1000; % number of total iterations
trans = 500; %number of initial points not to plot
P = zeros(tmax,1);

w = pmin:pstep:pmax;

%Redefine the map in terms of the parameter: 

fp = @(p,W) (1-W).*p + W.*(a-d-atan(lambda.*(p-c)))./b;

% initialise bifurcation diagram
fig3=figure(3);
fig2.Name = 'Bifurcation';
title('Bifurcation Diagram','FontSize',14,'interpreter','latex');
ylabel('$p_t$','FontSize',14,'interpreter','latex');
xlabel('$w$','FontSize',14,'interpreter','latex');

hold on
for k=1:length(w) % first loop iterates through parameter
    for i=2:tmax %same as before
        pp = P(i-1);
        P(i)=fp(pp,w(k));
    end
    % create a vector of w(k) in order to plot it:
    B = ones(1,tmax-trans)*w(k);
    scatter(B,P(trans+1:tmax),0.03,'k');
end


%% Parameter Basins for Periodic Cycles (Incomplete)
%{
cycle_num = 3; %number of periodic cycles 
cycles = 1:cycle_num;
max_dist = 1000;  % maximum distance after which the point is classified as 'divergent'

Px_min = 0; %minimum value for x-axis parameter
Px_max = 10;
Py_min = 0.1;
Py_max = 0.8;
step = 0.001;

tmax = 500; %number of iterations
trans = 200; %number of transient iterations
P = zeros(tmax,1);

init_min = 2; % initial condition min
init_max = 2; 
init_step = 0.1;
p0 = 4.8; %init_min:init_step:init_max;
P(1)=p0; 

lambda = Px_min:step:Px_max; %lambda is x-axis param.
w = Py_min:step:Py_max;

% Defining a colormap for up to 10 cycles
% white=non-convergent
% black = divergen
% red = steady-state
% 2-cycle = blue
cmap = [1 1 1; 0 0 0; 1 0 0; 0 0 1];

%Redefine the map in terms of the parameter: 
fp = @(p,L,W) (1-W).*p + W.*(a-d-atan(L.*(p-c)))./b;

[PX, PY]=meshgrid(lambda, w);
Basin = ones(numel(w),numel(lambda)); %Parameter basin matrix


for i=1:numel(PX)
    for i=2:tmax %same as before
        pp = P(i-1);
        P(i)=fp(pp,PX(i),PY(i));
        if P(i) - p0 > max_dist %break if divergent
            Basin(i) = -1; %-1 is color code for divergent
            break
        end
    end  
    C = numel(unique(P(trans:tmax))); 
    if C <= cycle_num
        Basin(i)= C;
    else
        Basin(i) = -2; % -2 is code for non-convergent
    end
end

fig4= figure(4);
[M,c]=contourf(lambda,w,Basin);
colormap(cmap);
%}
