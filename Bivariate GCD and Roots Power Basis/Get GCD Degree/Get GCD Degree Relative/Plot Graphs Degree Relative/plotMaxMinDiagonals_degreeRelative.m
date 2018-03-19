function [] = plotMaxMinDiagonals_degreeRelative(max_diags_R1, min_diags_R1, limits_k1, limits_k2, limits_t1, limits_t2)
%
% % Inputs
%
% max_diags_R1 : (Matrix)
%
% min_diags_R1 : (Matrix)
%
% limits_k1 : (Int Int)
%
% limits_k2 : (Int Int)
%
% limits_t1 : (Int Int)
%
% limits_t2 : (Int Int)

% Get upper and lower limits for k_{1} and k_{2}
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

% Get upper and lower limtis for t_{1} amd t_{2}
lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);
lowerLimit_t2 = limits_t2(1);
upperLimit_t2 = limits_t2(2);

% Get ratio
ratio = (log10(min_diags_R1) ./ log10(max_diags_R1));


% Plot
figure_name = sprintf([mfilename ' : ' 'MaxMin Row Diagonals' ]);
figure('name',figure_name)
hold on
x_vec = lowerLimit_k2 : 1 : upperLimit_k2;
y_vec = lowerLimit_k1 : 1 : upperLimit_k1;
[X, Y] = meshgrid(x_vec,y_vec);
surf(X,Y,(ratio));

hold off

end

