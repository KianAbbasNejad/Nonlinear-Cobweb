%% IterationMapFun
% This function calculates and plots the iteration loop for a discrete-time
% 1-dimensional map.
% (Kian Abbas Nejad)

% INPUTS
%   1. f - function handle: the map (function handle)
%   2. n - integer: the map's n-th iteration to plot symbolically 
%   3. x0 - number: initial condition
%   4. num - integer: number of iterations to plot from initial condition
%   4. name - string: name of the variable

% OPTIONAL INPUTS
%   5. ax - 4-vector: axis limits

function [] = IterationMapFun(f,n,x0,num,name,varargin)

% Calculating the n-th iteration of the map
fn = f; % creating new handle
for i=1:n-1  % recursive calculation of f^n(.)
    fn=@(x) f(fn(x));
end

figure('Name','IterationMap');
hold on
grid on
% Plotting the nth iteration of the map
fplot(fn,'r')
% add the 45-degree line
fplot(@(x) x,'b--')

% Plotting iterations from the initial condition
x1 = x0;
for i=1:num
    for k=1:n
        x1 =f(x1);
    end
    plot([x0 x0],[x0 x1],'k')
    plot([x0 x1],[x1 x1],'k')
    x0 = x1;
end

%% Graph Cleanup
ylab = ['$',name,'_{t+',num2str(n),'}','$'];
xlab = ['$',name,'_t','$'];
ylabel(ylab)
xlabel(xlab)
title(['Iteration of ',ylab])

% latex font and fontsize
a = gca;
set([a.Title a.XLabel a.YLabel],'Interpreter','Latex','FontSize',16);

% set axis limits
if nargin == 6
    axis(varargin{1})
else
    axis('fill')
end


end