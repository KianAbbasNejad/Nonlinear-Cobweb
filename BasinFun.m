%% BasinFun
% This code calculates the parameter basin for periodic cycles for a
% dynamical system.
% (Kian Abbas Nejad)
% INPUTS
%   1. f - function handle, with 3 inputs so that no.2 is for x-axis
%       parameter and no.3 is for y-axis parameter.
%   2. cycle_num - number of periodic cycles
%   3. xgrid - n-vector range of x-axis parameter
%   4. ygrid - n-vector range of y-axis parameter
%   5. init - initial condition
%   6. trans - number of transient iterations
%   7. itermax - number of total iterations
%   8. xname - string title for x-axis
%   9. yname - string title for y-axis
%
% OPTIONAL INPUTS as 1 str
%   10. 7. ax - 4-vector: axis limits

function [] = BasinFun(f,cycle_num,xgrid,ygrid,init,trans,itermax,xname,yname,varargin)
%% Algorithm Parameters
max_dist = 1000;  % maximum distance after which the point is classified as 'divergent'
tol = 3; %rounding tolerance '1e-tol' for number of periodic cycles

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

contour_levels= -1:cycle_num;

[PX, PY]=meshgrid(xgrid, ygrid);
Basin = -ones(numel(ygrid),numel(xgrid)); % Parameter basin matrix. -1 is non convergent

for i=1:numel(PX)
    P = [init; zeros(itermax,1)]; % time series
    for j=2:itermax 
        P(j)=f(P(j-1),PX(i),PY(i));
    end
    P = P(trans:itermax); % cutting transient iterations
    P = round(P,tol); % rounding using the tolerance level
    C = numel(unique(P));
    if (P(end) - P(1)) > max_dist % if divergent
        Basin(i) = 0; % 0 is color code for divergent
    elseif C <= cycle_num
        Basin(i)= C;
    end
end

% Plotting
figure('Name','Basin of Attraction');
hold on
[~,~]=contourf(xgrid,ygrid,Basin,contour_levels);
cmap = cmap(1:cycle_num+2,:);
colormap(cmap)

% Creating Legend
h = zeros(cycle_num+2,1);
for i=1:cycle_num+2
    h(i) = plot(NaN,NaN,'o','Color',cmap(i,:));
end
lgd = legend(h, labels(1:cycle_num+2),'Location','southoutside');
lgd.NumColumns = 4;

%% Graph Cleanup
title('Parameter Basin of Attraction');
xlab = ['$',xname,'$'];
ylab = ['$',yname,'$'];
ylabel(ylab);
xlabel(xlab);

% Latex font and fontsize
a = gca;
set([a.Title a.XLabel a.YLabel a.Legend],'Interpreter','Latex','FontSize',16);

% Set axis limits
if nargin == 10
    axis(varargin{1})
else
    axis('fill')
end

end