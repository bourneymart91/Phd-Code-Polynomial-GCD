function [residual,idx_col] = GetMinDistance(Sk)
% Given a Sylvester matrix remove each column in turn, and compute the
% distance between it and a vector which lies in the column space of the
% remaining columns.
% 
% Inputs.
%
% Sk : (Matrix) Sylvester subresultant matrix S_{k}(f,g)
%
% Outputs.
%
% residual : Residual obtained by removing the column ck. r = ck - Ax
%
% idx_col : (Int) Index of column removed where residual is minimal.



% Get the number of columns in the Sylvester matrix S_{k}(f,g)
nColumns = size(Sk,2);

% Initialise a vector to store residuals
vResiduals = zeros(nColumns,1);

for i = 1:1:nColumns
    
    % Get Ak, Sk with ck removed
    Ak = Sk;
    Ak(:,i) = [];
    
    % Get the column ck
    ck = Sk(:,i);
    
    % Get the solution x
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual
    vResiduals(i) = norm(ck - (Ak*x_ls));
    
end

% Get the minimal residual and the index of the corresponding column ck
[residual, idx_col] = min(abs(vResiduals));

%PlotMinimumResiduals(vResiduals);

end

function [] = PlotMinimumResiduals(vResiduals)


nEntries = length(vResiduals);
x_vec = 1 : 1 : nEntries;


figure_name = sprintf('Minimum Residuals');

figure('Name',figure_name)
hold on
plot(x_vec, log10(vResiduals))

xlabel('$i$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$ \log_{10} \left( r_{i} \right) $', 'Interpreter', 'latex', 'FontSize',20)
hold off

grid on
box on


end
