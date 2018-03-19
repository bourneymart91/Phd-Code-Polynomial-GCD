
close all; clc;




for mult = 1 : 1 : 10
    
    f = @(x,y) (x + 1) * (x + 2)^(mult) * (x - 3) * (x - 4) - y;
    g = @(x,y) -3* (x + 1) * (x + 2)^(mult) * (7*x + 11) - y ;
    figure()
    hold on
    grid on
    %f = @(x,y,z) (x^2 + y^2 + z^2 - 10); % sphere
    syms x y z
    interval_3D = [-2 2 -2 2 -2 2];
    interval_2D = [-3 0 -10 10];
    
    %fimplicit3(f,'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','red')
    %fimplicit3(g,'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','blue')
    
    %fimplicit3(f, interval_3D, 'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','red')
    %fimplicit3(g, interval_3D, 'EdgeColor', 'black', 'FaceAlpha', 0.8, 'FaceColor','blue')
    
    fimplicit(f, interval_2D);
    fimplicit(g, interval_2D);
    
    h1 = @(x,y) (x + 1);
    h2 = @(x,y) (x + 2);
    h3 = @(x,y) (x + 5);
    h4 = @(x,y) (x + 9);
    
    fimplicit(h1, interval_2D,'--g','LineWidth',2)
    fimplicit(h2, interval_2D,'--g','LineWidth',2)
    fimplicit(h3, interval_2D,'--g','LineWidth',2)
    fimplicit(h4, interval_2D,'--g','LineWidth',2)
    
    %fimplicit(h2, interval_2D, '--g','LineWidth',2);
    hold off
end