

close all; clc; 


for mult = 1:1:10

figure()
hold on 
grid on
%f = @(x,y,z) (x^2 + y^2 + z^2 - 10); % sphere
syms x y z


f = @(x,y,z) (x^2 + y^2 - 1) ^ mult  * (x + y + 0.2) - z;
g = @(x,y,z) (x^2 + y^2 - 1) ^ mult  * (y - 0.2) - z;

interval_3D = [-2 2 -2 2 -2 2];
interval_2D = [-2 2 -2 2];

%fimplicit3(f,'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','red')
%fimplicit3(g,'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','blue')

fimplicit3(f, interval_3D, 'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','red')
fimplicit3(g, interval_3D, 'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','blue')

h = @(x,y) (x.^2 + y.^2 - 1) ^ mult ;
fimplicit(h, interval_2D,'--g','LineWidth',2)

h2 = @(x,y) (x + 0.4);
fimplicit(h2, interval_2D, '--g','LineWidth',2);
hold off

end