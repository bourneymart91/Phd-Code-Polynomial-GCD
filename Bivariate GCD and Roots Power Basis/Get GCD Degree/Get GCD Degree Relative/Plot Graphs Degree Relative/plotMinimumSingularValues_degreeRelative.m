function [] = plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2)
%
% % Inputs
%
% mat_MinimumSingularValues : (Matrix)
%
% limits_k1 : (Int Int)
%
% limits_k2 : (Int Int)
%
% limits_t1 : (Int Int)
% 
% limits_t2 : (Int Int)

lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);
lowerLimit_t2 = limits_t2(1);
upperLimit_t2 = limits_t2(2);

figure_name = sprintf([ mfilename ' : ' 'Minimum Singular Values']);
figure('name', figure_name)

hold on
grid on

v_k1 = lowerLimit_k1 : 1 : upperLimit_k1;
v_k2 = lowerLimit_k2 : 1 : upperLimit_k2;

[X,Y] = meshgrid(v_k1,v_k2);

mesh(X, Y, log10(mat_MinimumSingularValues)');


z_min = min(min(log10(mat_MinimumSingularValues)));
z_max = max(max(log10(mat_MinimumSingularValues)));

points_x_lower = [lowerLimit_t1 lowerLimit_t1 lowerLimit_t1 lowerLimit_t1 ];
points_x_upper = [upperLimit_t1 upperLimit_t1 upperLimit_t1 upperLimit_t1];
points_y = [lowerLimit_k2 upperLimit_k2 upperLimit_k2 lowerLimit_k2 ];
points_z = [z_max z_max z_min z_min];

lower_bound_x_patch = patch(points_x_lower, points_y, points_z,'red', 'facealpha', 0.3);
upper_bound_x_patch = patch(points_x_upper, points_y, points_z,'red', 'facealpha', 0.3);

points_y_lower = [lowerLimit_t2 lowerLimit_t2 lowerLimit_t2 lowerLimit_t2];
points_y_upper = [upperLimit_t2 upperLimit_t2 upperLimit_t2 upperLimit_t2];
points_x = [lowerLimit_k1, upperLimit_k1  upperLimit_k1 lowerLimit_k1];

lower_bound_y_patch = patch(points_x, points_y_lower, points_z,'blue', 'facealpha', 0.3);
upper_bound_y_patch = patch(points_x, points_y_upper, points_z,'blue', 'facealpha', 0.3);



xlabel('k_{2}')
ylabel('k_{1}')
zlabel('log_{10} \sigma_{k1,k2}')

hold off

end