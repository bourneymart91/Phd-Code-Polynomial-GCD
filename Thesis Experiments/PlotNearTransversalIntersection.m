function [] = PlotNearTransversalIntersection()

% Get two curves whose intersection is transverse
figure('name','Plot Transverse Intersection')
fplot(@(x) 3*x +1);
hold on
grid on
fplot(@(x) 3.2*x +1);
hold off

end