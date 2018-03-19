function [] = o_intersections

% Define an Implicit surface
syms x y z
imp_surface = x^2 + y^2 + z^2;



bez_surf_x = ...
    [...
    1 2 3 4
    2 3 4 0
    3 4 0 0
    4 0 0 0
    ];
vBez_surf_x = GetAsVector(bez_surf_x)

bez_surf_y = ...
    [...
    1 2 3 4
    1 2 3 0
    1 2 0 0
    1 0 0 0
    ];
vBez_surf_y = GetAsVector(bez_surf_y)

bez_surf_z = ...
    [
    1 2 3 4
    1 2 3 0
    1 2 0 0
    1 0 0 0
    ];
vBez_surf_z = GetAsVector(bez_surf_z)

figure()
figure()
hold on
ezimplot3('x^2 + y^2 + z^2 - 5 ')
scatter3(vBez_surf_x,vBez_surf_y,vBez_surf_z)
hold off



end

