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
n = 1; %number of iterations
p0 = -6; %initial condition
p1 = p0; %vector for iteration loop
num = 5; % number of iterations for plotting
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

%% Bifurcation Diagram
%{
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
%}

%% Color Map for Parameter Basins
% Defining a colormap for up to 10 cycles
% non-convergent(-1) = white    divergent(0) = black
% steady-state  = red           2-cycle = blue
% 3-cycle = yellow              4-cycle = orange
% 5-cycle = gray                6-cycle = teal
% 7-cycle = purple              8-cycle = olive
% 9-cycle = green               10-cycle = brown
cmap = [1 1 1; 0 0 0
    1 0 0; 0 0 1
    1 1 0; 1 0.5 0
    0.5 0.5 0.5; 0 0.5 0.5
    0.5 0 0.5 ; 0.5 0.5 0
    0 0.5 0; 0.5 0.27 0.1];
% labels for legend:
labels = {'Non-convergent','divergent','steady-state','2-cycle','4-cycle',...
    '5-cycle','6-cycle','7-cycle','8-cycle','9-cycle','10-cycle'};
%% Parameter Basins for Periodic Cycles 

cycle_num = 7; %number of periodic cycles 
contour_levels= -1:cycle_num; 
max_dist = 1000;  % maximum distance after which the point is classified as 'divergent'
tol = 3; %rounding tolerance '1e-tol' for number of periodic cycles
p0 = 6; %initial condition
tmax = 1000; %number of iterations
trans = 700; %number of transient iterations

Px_min = 0; %minimum value for x-axis parameter (here lambda)
Px_max = 10;
Py_min = 0.1;
Py_max = 0.8;
step = 0.002;

lambda = Px_min:step:Px_max; %lambda is x-axis param.
w = Py_min:step:Py_max;

%Redefine the function in terms of the parameter: 
Fp = @(p,L,W) (1-W).*p + W.*(a-d-atan(L.*(p-c)))./b;

[PX, PY]=meshgrid(lambda, w);
Basin = -ones(numel(w),numel(lambda)); %Parameter basin matrix. -1 is non convergent

for i=1:numel(PX)
    P = [p0; zeros(tmax,1)];
    for j=2:tmax %same as before
        pp = P(j-1);
        P(j)=Fp(pp,PX(i),PY(i));
    end  
    P=round(P,tol); %rounding using the tolerance level
    C = numel(unique(P(trans:tmax))); 
    if (P(end) - P(1)) > max_dist % if divergent
        Basin(i) = 0; %0 is color code for divergent
    elseif C <= cycle_num
        Basin(i)= C;
    end 
end
fig4= figure(4);
hold on
[M,c]=contourf(lambda,w,Basin,contour_levels);
cmap = cmap(1:cycle_num+2,:);
colormap(cmap)
fig2.Name = 'Basin_of_Attraction';
title('Parameter Basin of Attraction','FontSize',14,'interpreter','latex');
ylabel('$w$','FontSize',14,'interpreter','latex');
xlabel('$\lambda$','FontSize',14,'interpreter','latex');

% Creating Legend
h = zeros(cycle_num+2,1);
for i=1:cycle_num+2
    h(i) = plot(NaN,NaN,'o','Color',cmap(i,:));
end
lgd = legend(h, labels(1:cycle_num+2),'Location','southoutside','interpreter','latex');
lgd.NumColumns = 4;
