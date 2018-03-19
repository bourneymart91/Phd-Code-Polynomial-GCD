function [] = PlotExplicitSurface(fxy)
% Given the power basis coefficients, plot the surface z = f(x,y)


% % Evaluate the power surface over the interval

lwr_bound_x = -2.5;
upr_bound_x = 0;

lwr_bound_y = -1.2;
upr_bound_y = -0.8;

inc =0.1;

x_val_vec = lwr_bound_x : inc : upr_bound_x;
y_val_vec = lwr_bound_y : inc : upr_bound_y;

[p,q] = meshgrid(x_val_vec, y_val_vec);

for i = 1:1:length(x_val_vec)
    for j = 1:1:length(y_val_vec)
        z(i,j) = Evaluate_PowerPoly_Bivariate(x_val_vec(i),y_val_vec(j),...
                        fxy);
    end
end
figure('name','Polynomial Plot')
surf(p,q,z')
hold off



end