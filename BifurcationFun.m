%% BifurcationFun 
% Calculates and plots the bifurcation diagram of a map
% (Kian Abbas Nejad)
% INPUTS
%   1. f - function handle: with 2 inputs, second one being the parameter of
%       the dynamic system 
%   2. prange - n-vector of parameter range
%   3. trans - number of transient iterations
%   4. itermax - number of total iterations 
%   5. xname - string: name of dynamic system
%   6. pname - string: name of parameter

% OPTIONAL INPUTS
%   7. ax - 4-vector: axis limits

function [] = BifurcationFun(f,prange,trans,itermax,xname,pname,varargin)


P = zeros(itermax,1);

% Initialise bifurcation diagram
figure('Name','Bifurcation');
hold on

for k=1:length(prange) % first loop iterates through parameter
    for i=2:itermax %same as before
        pp = P(i-1);
        P(i)=f(pp,prange(k));
    end
    % create a vector of w(k) in order to plot it:
    B = ones(1,itermax-trans)*prange(k);
    scatter(B,P(trans+1:itermax),0.03,'k');
end


%% Graph Cleanup
title('Bifurcation Diagram');

xlab = ['$',pname,'$'];
ylab = ['$',xname,'_t','$'];
ylabel(ylab);
xlabel(xlab);

% latex font and fontsize
a = gca;
set([a.Title a.XLabel a.YLabel],'Interpreter','Latex','FontSize',16);

% set axis limits
if nargin == 7
    axis(varargin{1})
else
    axis('fill')
end

end