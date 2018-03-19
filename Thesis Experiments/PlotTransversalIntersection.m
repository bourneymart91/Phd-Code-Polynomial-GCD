function [] = PlotTransversalIntersection()

% Get two curves whose intersection is transverse
figure('name','Plot Transverse Intersection')
fplot(@(x) 2*x +2);
hold on
fplot(@(x) -2*x +2);
hold off

end