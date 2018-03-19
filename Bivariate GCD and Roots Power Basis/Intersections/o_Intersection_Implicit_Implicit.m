function [] = o_Intersection_Implicit_Implicit()
%
% Compute the intersections between two implicitly defined curves.
%
% o_Intersection_Implicit_Implicit()
ex_num = '2';

switch ex_num
    case '1'
        
        p = ...
            [
            0 0 1
            0 0 0
            -1 0 0
            1 0 0
            ];
        
        q = ...
            [
            0 0 1
            -1 0 0
            2 0 0
            -1 0 0
            ];
        
        
        m = 3;
        n = 3;
        f = polynomial_subtraction(p,q);
        
    case '2'
        p = ...
            [
            0 0 1 -1 1
            0 0 0 0 0
            0 -2 0 0 0
            0 0 0 0 0
            1 0 0 0 0
            ];
        q = ...
            [
            0 1 
            0 0
            -2 0
            ];
        m = 4;
        n = 2;
        
        f = polynomial_subtraction(p,q);
        
        
            
end


% %
% %

problem_type = 'Intersection';
ex_num = 'Intersection';
emin = 0;
emax = 0;
mean_method = 'None';
bool_alpha_theta = 'n';
low_rank_approx_method = 'None';
SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method)

% %
% %


o_roots_mymethod(f,m)





[fxy, gxy, dxy, uxy, vxy, t, t1, t2] = o_gcd_mymethod(p,q,3,3,[0,min(m,n)]);

end

function h = polynomial_subtraction(fxy, gxy)

[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

max_x_power = max(m1,n1);
max_y_power = max(m2,n2);

temp = zeros(max_x_power+1,max_y_power+1);

fxy = temp;
fxy(1:m1+1,1:m2+1) = fxy;
gxy = temp;
gxy(1:n1+1,1:n2+1) = gxy;

h = fxy - gxy;

end