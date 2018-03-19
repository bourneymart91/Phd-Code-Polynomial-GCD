function [lambda, mu, rho, theta] =  OptimalAlphaBetaGammaTheta_3Polys_3Eqns(...
    vF_max1, vF_min1, ...
    vF_max2, vF_min2, ...
    vG_max1, vG_min1, ...
    vG_max2, vG_min2, ...
    vH_max1, vH_min1, ...
    vH_max2, vH_min2 ...
    )
% OptimalAlphaBetaTheta(F_max,F_min,G_max,G_min)
%
% This function computes the optimal values alpha and theta for the
% preprocessing operations on the kth Sylvester subresultant matrix
% S_{k}(f,g)Q_{k} = D^{-1}_{k}T_{k}Q_{k}
%
%
% % Inputs
%
%
% vF_max :  (Vector) A vector of length m+1, such that F_max(i) stores the
%          element of maximum magnitude of S(f,g)Q that contains the
%          coefficient a(i) of f, i=1,...,m+1.
%
% vF_min :  (Vector) A vector of length m+1, such that F_min(i) stores the
%          element of minimum magnitude of S(f,g)Q that contains the
%          coefficient a(i) of f, i=1,...,m+1.
%
% vF_max2 : (Vector)
%
%
% vF_min2 : (Vector)
%
%
% vG_max : (Vector) A vector of length n+1, such that G_max(i) stores the
%         element of maximum magnitude of S(f,g)Q that contains the
%         coefficient b(i) of g, i=1,...,n+1.
%
% vG_min : (Vector) A vector of length n+1, such that G_min(i) stores the
%         element of minimum magnitude of S(f,g)Q that contains the
%         coefficient b(i) of g, i=1,...,n+1.
%
% vH_max : (Vector)
%
% vH_min : (Vector)
%
% % Outputs
%
% alpha :(Float)
%
% beta : (Float)
%
% theta : (Float)
%



% Define vector f

global SETTINGS
method = SETTINGS.SCALING_METHOD;

% get degree of polynomials f and g
m = GetDegree(vF_max1);
n = GetDegree(vG_max1);
o = GetDegree(vH_max1);


switch method
    case 'lambda_mu_rho'
        
        f = [1, -1, 0, 0, 0, 0];
        
        % Assemble the four submatrices.
        A1 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'   -ones(m+1,1)     zeros(m+1,1)  zeros(m+1,1) ];
        A2 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'   zeros(n+1,1)     -ones(n+1,1)  zeros(n+1,1) ];
        A3 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'   -ones(m+1,1)     zeros(m+1,1)  zeros(m+1,1) ];
        A4 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'   zeros(o+1,1)     zeros(o+1,1)  -ones(o+1,1) ];
        A5 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'   zeros(o+1,1)     zeros(o+1,1)  -ones(o+1,1) ];
        A6 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'   zeros(n+1,1)     -ones(n+1,1)  zeros(n+1,1) ];
        
        A7 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'    ones(m+1,1)      zeros(m+1,1)  zeros(m+1,1) ];
        A8 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'    zeros(n+1,1)     ones(n+1,1)   zeros(n+1,1) ];
        A9 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'    ones(m+1,1)      zeros(m+1,1)  zeros(m+1,1) ];
        A10 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'    zeros(o+1,1)     zeros(o+1,1)   ones(o+1,1) ];
        A11 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'    zeros(o+1,1)     zeros(o+1,1)   ones(o+1,1) ];
        A12 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'    zeros(n+1,1)     ones(n+1,1)   zeros(n+1,1) ];
     
    case 'lambda_mu'
        f = [1, -1, 0, 0, 0];
        
        % Assemble the four submatrices.
        A1 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'   -ones(m+1,1)     zeros(m+1,1)   ];
        A2 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'   zeros(n+1,1)     -ones(n+1,1)   ];
        A3 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'   -ones(m+1,1)     zeros(m+1,1)   ];
        A4 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'   zeros(o+1,1)     zeros(o+1,1)   ];
        A5 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'   zeros(o+1,1)     zeros(o+1,1)   ];
        A6 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'   zeros(n+1,1)     -ones(n+1,1)   ];
        
        A7 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'    ones(m+1,1)      zeros(m+1,1)   ];
        A8 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'    zeros(n+1,1)     ones(n+1,1)    ];
        A9 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'    ones(m+1,1)      zeros(m+1,1)   ];
        A10 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'    zeros(o+1,1)     zeros(o+1,1)  ];
        A11 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'    zeros(o+1,1)     zeros(o+1,1)  ];
        A12 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'    zeros(n+1,1)     ones(n+1,1)   ];
     
    case 'mu_rho'
        
        f = [1, -1, 0, 0, 0];
        
        % Assemble the four submatrices.
        A1 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'        zeros(m+1,1)  zeros(m+1,1) ];
        A2 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'        -ones(n+1,1)  zeros(n+1,1) ];
        A3 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'        zeros(m+1,1)  zeros(m+1,1) ];
        A4 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'        zeros(o+1,1)  -ones(o+1,1) ];
        A5 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'        zeros(o+1,1)  -ones(o+1,1) ];
        A6 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'        -ones(n+1,1)  zeros(n+1,1) ];
        
        A7 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'          zeros(m+1,1)  zeros(m+1,1) ];
        A8 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'          ones(n+1,1)   zeros(n+1,1) ];
        A9 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'          zeros(m+1,1)  zeros(m+1,1) ];
        A10 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'         zeros(o+1,1)   ones(o+1,1) ];
        A11 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'         zeros(o+1,1)   ones(o+1,1) ];
        A12 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'         ones(n+1,1)   zeros(n+1,1) ];
    
    case 'lambda_rho'
    
        f = [1, -1, 0, 0, 0];
        
        % Assemble the four submatrices.
        A1 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'   -ones(m+1,1)       zeros(m+1,1) ];
        A2 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'   zeros(n+1,1)       zeros(n+1,1) ];
        A3 = [ ones(m + 1,  1)   zeros(m + 1, 1) -(0:1:m)'   -ones(m+1,1)       zeros(m+1,1) ];
        A4 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'   zeros(o+1,1)       -ones(o+1,1) ];
        A5 = [ ones(o + 1,  1)   zeros(o + 1, 1) -(0:1:o)'   zeros(o+1,1)       -ones(o+1,1) ];
        A6 = [ ones(n + 1,  1)   zeros(n + 1, 1) -(0:1:n)'   zeros(n+1,1)       zeros(n+1,1) ];
        
        A7 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'    ones(m+1,1)        zeros(m+1,1) ];
        A8 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'    zeros(n+1,1)        zeros(n+1,1) ];
        A9 = [ zeros(m + 1, 1)  -ones(m + 1, 1)  (0:1:m)'    ones(m+1,1)         zeros(m+1,1) ];
        A10 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'    zeros(o+1,1)        ones(o+1,1) ];
        A11 = [ zeros(o + 1, 1)  -ones(o + 1, 1)  (0:1:o)'    zeros(o+1,1)        ones(o+1,1) ];
        A12 = [ zeros(n + 1, 1)  -ones(n + 1, 1)  (0:1:n)'    zeros(n+1,1)       zeros(n+1,1) ];
     
        
        
    otherwise
        error('err')
end



A = -[A1; A2; A3; A4; A5; A6; A7; A8; A9; A10; A11 ; A12];

% Construct the vector b.
b= -[
    log10(abs(vF_max1));
    log10(abs(vG_max1));
    log10(abs(vF_max2));
    log10(abs(vH_max1));
    log10(abs(vH_max2));
    log10(abs(vG_max2));
    -log10(abs(vF_min1));
    -log10(abs(vG_min1));
    -log10(abs(vF_min2));
    -log10(abs(vH_min1));
    -log10(abs(vH_min2));
    -log10(abs(vG_min2));
    ];

% Solve the linear programming problem and extract alpha and theta
% from the solution vector x.

warning('off')
x = linprog(f, A, b);
warning('on')

switch method
    case 'lambda_mu_rho'

        theta = real(10^x(3));
        lambda = real(10^x(4));
        mu = real(10^x(5));
        rho = real(10^x(6));

    case 'lambda_mu'    
        
        theta = real(10^x(3));
        lambda = real(10^x(4));
        mu = real(10^x(5));
        rho = 1;
        
    case 'mu_rho'
        
        theta = real(10^x(3));
        lambda = 1;
        mu = real(10^x(4));
        rho = real(10^x(5));

    case 'lambda_rho'
        
        theta = real(10^x(3));
        lambda = 1;
        mu = real(10^x(4));
        rho = real(10^x(5));
        
        
        
end




end