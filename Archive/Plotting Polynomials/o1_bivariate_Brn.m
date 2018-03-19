function []= o1_bivariate_Brn()


example_style = input('Example Style From Separable roots (r) or coefficients (c):  ','s')

switch example_style
    case 'r'
        
        % Input the example Number
        ex_num = input('Example Number: ','s');
        
        % Get the roots of polynomial f
        [f_x_roots,f_y_roots] = Examples_Bivariate_Separable(ex_num);
        
        % Get the coefficients of the polynomials f(x) and f(y)
        f_x_poly_bi = BuildPoly_Brn(f_x_roots);
        f_y_poly_bi = BuildPoly_Brn(f_y_roots);
        
        % Get the degrees m1 and m2 of polynomials f wrt x and f wrt y,
        % respectively.
        m1 = length(f_x_poly_bi)-1;
        m2 = length(f_y_poly_bi)-1;
        
        % Get matrix of coefficients of f(x,y) including binomial coefficients.
        fxy_matrix_bi_Brn = f_x_poly_bi * f_y_poly_bi';
        
        % Remove binomial coefficients corresponding to x, from the rows of fxy
        for i = 0:1:m1
            fxy_matrix_Brn(i+1,:) = fxy_matrix_bi_Brn(i+1,:) ./ nchoosek(m1,i);
        end
        % Remove binomial coefficients corresponding to y, from the columns of fxy
        for j = 0:1:m2
            fxy_matrix_Brn(:,j+1) = fxy_matrix_Brn(:,j+1) ./ nchoosek(m2,j)  ;
        end
        
    case 'c'
        ex_num = input('Enter Example Number :  ','s')
        
        switch ex_num
            case '1'
                %fxy_matrix_Brn = [2 2 2; 2 -1 3; 2 2 2];
                fxy_matrix_Brn = ...
                    [  
                        2    2   2   2;
                        2   -3  -3   2;
                        2   -3  -3   2;
                        2    2   2   2;
                    ];
                
            case '2'
                
                fxy_matrix_Brn = ...
                    [
                        -1   1   1   1;
                         1   1   1   1;
                         1   1   2   1;
                         1   1   1   1;
                    ]
        end
        
        % Get the degrees m1 and m2 of polynomials f wrt x and f wrt y,
        % respectively.
        rows = size(fxy_matrix_Brn,1)  ;
        cols = size(fxy_matrix_Brn,2)  ;
        m1 = rows -1;
        m2 = cols -1;
        
        
end


% % % Given that we have obtained coefficients for fxy in the power and
% % bernstein basis, we now add the same noise to both sets of coefficients.
%
noise = input('Do you wish to add noise? (Y)es or (N)o :  ','s')
switch noise
    case 'N'
        % Set upper and lower noise threshold
        emin = 1e-8;
        emax = 1e-2;
        
        % add noise to the coefficients in Bernstein basis
        fxy_matrix_Brn = Noise(fxy_matrix_Brn,emin,emax);
        
        
    case 'Y'
end



%% Plot 1
% % Plot the control points of the surface in bernstein basis
rows = size(fxy_matrix_Brn,1)
cols = size(fxy_matrix_Brn,2)

fxy_Coordinates_CP = [];
for i = 0:1:rows-1
    for j = 0:1:cols-1
        % get x coordinate of control point
        x_pos = i/(rows-1);
        % Get y coordinate of control point
        y_pos = j/(cols-1);
        % Get z coordinate of control point
        z_pos = fxy_matrix_Brn(i+1,j+1);
        
        fxy_coord = [x_pos y_pos z_pos];
        
        % Add the coordinates to the dataset
        fxy_Coordinates_CP = [fxy_Coordinates_CP; fxy_coord];
    end
end

% Plot the control points in a three dimensional scatter
figure(1)
hold on
scatter3(fxy_Coordinates_CP(:,1), fxy_Coordinates_CP(:,2), fxy_Coordinates_CP(:,3),'Filled')
hold off

%% Plot 2
% % Plot the control points of the surface as a surface.
steps_x = 0:1/(rows-1):1;
steps_y = 0:1/(cols-1):1;

[XI_p1,YI_p1] = meshgrid(steps_x, steps_y);

% now interpolate - find z values for these grid points
ZI_p1 = griddata(fxy_Coordinates_CP(:,1),fxy_Coordinates_CP(:,2),fxy_Coordinates_CP(:,3),XI_p1, YI_p1);


% Plot the control polygon
figure(2)
g1 = mesh(XI_p1,YI_p1,ZI_p1,'edgecolor','black');
hidden off;
hold off

%%
% Split the x variable into m1 chunks
fxy_coordinates_Surface = []
for i = 0:1:m1
    x_inc = i/m1;
    for j = 0:1:m2
        y_inc = j/m2;
        z_Brn_exact = Evaluate_BernsteinPoly(x_inc,y_inc,...
            fxy_matrix_Brn);
        fxy_coordinates_Surface = [fxy_coordinates_Surface; x_inc y_inc z_Brn_exact]
        
    end
end

%% Plot 3
% % Plot the Bézier surface by evaluating at a series of points

fxy_coordinates_Surface = [];

% % Evaluate the Bernstein surface over the interval
for i = 0:0.1:m1
    x_inc = i/m1;
    for j = 0:0.1:m2
        y_inc = j/m2;
        
        z_Brn_exact = Evaluate_BernsteinPoly(x_inc,y_inc,...
            fxy_matrix_Brn);
        
        fxy_coordinates_Surface = [fxy_coordinates_Surface; x_inc y_inc z_Brn_exact];
        
    end
end


% % plot the exact Bernstein data and the noisy bernstein data
figure(3)
title('Bernstein Surface')
hold on

% Plot exact data
scatter3(fxy_coordinates_Surface(:,1),fxy_coordinates_Surface(:,2),fxy_coordinates_Surface(:,3),'Filled')

hold off

%% Plot 4
% % Plot the bezier surface from the points above
% Get x y and z components of the exact Bernstein data
x_Brn_exact = fxy_coordinates_Surface(:,1);
y_Brn_exact = fxy_coordinates_Surface(:,2);
z_Brn_exact = fxy_coordinates_Surface(:,3);

% Get step size
steps = 0:.05:1;

[XI,YI] = meshgrid(steps, steps);
% now interpolate - find z values for these grid points
ZI_plot4 = griddata(x_Brn_exact,y_Brn_exact,z_Brn_exact,XI, YI);

% % Plot the exact and noisy bernstein surface
figure(4)
title('Exact Bernstein Surface')
hold on
g1 = surf(XI,YI,ZI_plot4);
g3 = mesh(XI_p1,YI_p1,ZI_p1,'edgecolor','black','facecolor','blue','MarkerFaceColor','black');
alpha(g3,0.05)
hold off



end



