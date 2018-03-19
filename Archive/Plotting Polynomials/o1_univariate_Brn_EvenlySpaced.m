function [] = o1_univariate_Brn_EvenlySpaced()
% Given a univariate polynomials roots, build the polynomial, plot the control
% points, draw the curve. Bernstein Basis.

example_style = input('Example Style : From Roots (r) or From Coefficients (c):  ','s');

switch example_style
    case 'r'
        %% Build a Polynomial whose coefficients give the y value of the Bernstein Polynomial
        
        % Get the set of roots
        [roots, ~] = Univariate_Roots_Examples(1)
        
        % Get coefficients of f including binomial coefficients.
        f_bi = Build_poly(roots);
        
        % Get degree of polynomial f
        m = size(f_bi,1)-1;
        
        % Build matrix of binomial coefficients.
        Bi_m = zeros(m+1,1);
        for i=0:1:m
            Bi_m(i+1) = nchoosek(m,i);
        end
        
        % Get fx by removing binomial coefficients from f_bi
        fx = f_bi./Bi_m
        
        % Get control points
        CP = GetControlPoints(0,1,fx);
        
    case 'c'
        % % Build a polynomial whose coefficients give the y value of the Bernstein Polynomial
        % Provide a polynomial defined by coefficients in the Bernstein
        % Basis
        fx = [1 2 1 3]
        CP = GetControlPoints(0,1,fx);
        
        
end

%% Plot Figure 1 - Produce a plot of the Curve defined by the control points


% produce a set of points to evaluate fx
for i = 1:1:100
    y_data(i) = Bernstein_Eval(fx,i/100);
    x_data(i) = i/100;
end

% Plot the control points and the data.
figure(1)
hold on
plot(x_data,y_data,'-')
plot(CP(:,1),CP(:,2),'-s')
hold off

end


function Pk = GetControlPoints(a,b,f)
% Get set of control points for the function f on the interval [a,b]
% a :- Interval Lower Limit.
% b :- Interval Upper Limit.
% f :- Coefficients of Polynomial
%
% Pk :- Set of control points.

% Get degree of polynomial f
m = length(f)-1;

% Initialise empty vector of control points.
Pk = [];

% for each control point, assign value.
for i = 0:1:m
    Pk = [Pk; a+(i/m).*(b-a)    f(i+1)];
end



end

function [f_bi]=Build_poly(A)
% Obtain polynomial coefficients in the Scaled Bernstein Basis, given a set of roots and multiplicities.
% This function implements the convolution operation on the polynomial
% defined by the matrix A, which has two columns, that is, this function
% computes the coefficients of the polynomial defined by A. These are
% the coefficients in the scaled Bernstein basis form of the polynomial.
% Coefficients of the form ai(m choose i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of distinct roots of the polynomial.
r = size(A,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f_bi = 1;
for k = 1:1:r
    w = B_conv(A(k,1),A(k,2));
    f_bi = conv(f_bi,w) ;
end
f_bi = f_bi';
end


function [t]=B_conv(root,mult)
% This function convolves the vector [-r 1-r] with itself m times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs:


% r :   Root

% m :   Multiplicity of root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs:


% t :   Vector which stores the result from this convolution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the
% scaled Bernstein basis.


if mult==1
    t=[-root,1-root];
else
    
    q=[-root,1-root];
    t=[-root,1-root];
    for k=2:1:mult
        t=conv(t,q);
    end
end
end





function [xx] =  Bernstein_Eval(f,c)
% Evaluate function f at point c.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% f :   Polynomial f

% c :   Point at which we evaluate polynomial f.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = c;
m = length(f)-1;
x = zeros(length(f),1);

for i = 0:length(f)-1
    x(i+1) = f(i+1) .* nchoosek(m,i) .* (1-y)^(m-i) .* y^i ;
end

xx = sum(x);

end
