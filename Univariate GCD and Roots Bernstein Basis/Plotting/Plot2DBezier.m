function [] = Plot2DBezier(cp_f_arr)
% Plot the curves C_{1},C_{2},... whose control points are given in an
% array cp_f_arr.
%
% Inputs.
% 
% cp_f_arr : (Array of Matrices) Cell array, where each cell cp_f_arr{i} contains the control
% points of curve C_{i}
%
%
% Plot2DBezier({[1 2.5 6.5 7; 2 5 6 3]})


% matrix of control points must contain 2 rows and m+1 columns


% Plot the Graph
figure('name','Bezier Plot')
hold on

% Get the number of sets of control points
[~,nCurves] = size(cp_f_arr);

for curve_number = 1 : 1 : nCurves
    
    % Get the current set of control points
    controlPoints_fx = cp_f_arr{curve_number};
    
    % Given a set of control points, plot the bezier curve.
    degree_fx = size(controlPoints_fx, 2) - 1;
    
    t = linspace(0, 1, 101);
    
    pts_f = 0;
    for i = 0 : 1 : degree_fx
        pts_f = pts_f + kron( nchoosek(degree_fx,i).*((1-t).^(degree_fx - i)) .* (t.^(i)) ,controlPoints_fx(:,i+1));
    end

    % For every column (control point) of f
    for i = 1:1:degree_fx+1
        str = sprintf('c_%i',i);
        placelabel(controlPoints_fx(:,i),str);
    end
    
    plot_name = sprintf('Curve C_%i',curve_number);
    plot(pts_f(1,:),pts_f(2,:),'DisplayName',plot_name)
    
    
       
end

grid on

legend('show')
hold off

end

function placelabel(pt,str)
x = pt(1);
y = pt(2);
h = line(x,y);
h.Marker = '.';
h = text(x,y,str);
h.HorizontalAlignment = 'center';
h.VerticalAlignment = 'bottom';
end