function [] = o_intersection_Bernstein_Bernstein(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
% Given two Bernstein polynomials, calculate the points of intersection
%
% Inputs.
%
% ex_num : Example Number
%
% emin : Lower noise level
%
% emax : Upper noise level
%
% mean_method :
%
% bool_alpha_theta :
%
% low_rank_approx_method :
%
% apf_method
%
%
% % Example
% 
% o_intersection_Bernstein_Bernstein('1',1e-10,1e-12,'Geometric Mean Matlab Method','y','Standard STLN','None')

% Set global variables
SetGlobalVariables(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method, apf_method)

% Get the two example Curves
[C1,C2] = Examples_Explicit_Bezier_Intersection(ex_num);

% Plot the curves
Plot2DBezier({C1,C2});

% get f(x) = y and g(x) = y
fx = C1(2,:)';
gx = C2(2,:)';

% get f(y) = x and g(y) = x
fy = C1(1,:)';
gy = C2(1,:)';

% Add Noise

% Get the implicit equation, and find the roots
hx = fx - gx;

% % Get Roots by Matlab Method
roots_matlab = o_roots_matlab(hx)

%
roots_mymethod = o_roots_mymethod(hx)

%%
% Given the roots get the points of intersection.
% Evaluate the Bézier curve at the points of intersection to find y
% coordinate.
% % 
% Given roots are calculated, plug in roots to one of the original 
% parametric equations
[num_rts,~] = size(roots_matlab);

% for each root
for i = 1:1:num_rts
    r_i = roots_matlab(i);
    if isreal(r_i)
        
        
        x1 = Bernstein_Evaluate(fy,r_i);
        y1 = Bernstein_Evaluate(fx,r_i);
        x2 = Bernstein_Evaluate(gy,r_i);
        y2 = Bernstein_Evaluate(gx,r_i);
        
        display([x1 y1])
        display([x2 y2])
    %fprintf('The root %2.4e gives intersection point %2.3f , %2.3f \n',r_i,x,y)
   
    end
end
end
