function [] = plot_fxy(fxy_matrix_Pwr)
% Given the power basis coefficients, plot the polynomial fxy


% % Evaluate the power surface over the interval

lwr_bound_x = 0;
upr_bound_x = 3.4;

lwr_bound_y = -3;
upr_bound_y = 1;

inc =0.1;

x_val_vec = lwr_bound_x : inc : upr_bound_x;
y_val_vec = lwr_bound_y : inc : upr_bound_y;

[p,q] = meshgrid(x_val_vec, y_val_vec);
pairs = [p(:) q(:)];

% initialise a vector to store z values (exact)
r = zeros(size(pairs,1),1);

% For each pair of (x,y) values
for i = 1:1:size(pairs,1)
    
    % Get the x,y value
    pair = pairs(i,:);
    x_val = pair(1);
    y_val = pair(2);
    
    % Get the z value
    z_val           =   Evaluate_PowerPoly_Bivariate(x_val,y_val,...
                        fxy_matrix_Pwr);
         
                    
    % Add the z value to the vector r
    r(i) = z_val;

                    
end


data_pwr = [pairs r];

% % Plot the exact and noisy power basis surface
figure(3)
title('Exact and Noisy Power Surface')
hold on
scatter3(data_pwr(:,1),data_pwr(:,2),data_pwr(:,3),'filled')
hold off

% Get the x y and z coordinates of the exact power coefficients
x_vec = data_pwr(:,1);
y_vec = data_pwr(:,2);
z_vec = data_pwr(:,3);


% Get step size
steps_x = lwr_bound_x : inc : upr_bound_x;
steps_y = lwr_bound_y : inc : upr_bound_y;

[XI,YI] = meshgrid(steps_x, steps_y);
% now interpolate - find z values for these grid points
ZI = griddata(x_vec,y_vec,z_vec,XI, YI);

figure(1)
title('Surface Plotting')
hold on
xlabel('x')
ylabel('y')
zlabel('z')
g1 = surf(XI,YI,ZI);
alpha(g1,0.5);
hold off


end