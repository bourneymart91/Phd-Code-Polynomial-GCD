function [lambda, mu, rho, theta] = GetOptimalAlphaBetaAndTheta_3Polys(fx, gx, hx)
% Compute the optimal values of \alpha, \beta and \theta
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% hx : (Vector) Coefficients of polynomial h(x)
%
% % Outputs
%
% alpha : (Float)
%
% beta : (Float) 
%
% theta : (Float)

% Ensure that f(x) and g(x) are column vectors
if (size(fx,2) >1  || size(gx,2)>2 || size(gx,2) > 1)
    error('f(x), g(x) and h(x) must be column vectors')
end



% Get degree of polynomial f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

scaling_method = 'lambda_rho';

switch scaling_method

    case 'lambda_mu_rho'
    
        f = [1 -1 0  0 0 0 ];
        % Build the first partiion
        A1 = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)'  -ones(m + 1, 1)     zeros(m + 1, 1)         zeros(m + 1, 1) ];
        A2 = [ones(n+1,1) zeros(n+1,1) -1.*(0:1:n)'  zeros(n + 1, 1)     -ones(n + 1, 1)         zeros(n + 1, 1)];
        A3 = [ones(o+1,1) zeros(o+1,1) -1.*(0:1:o)'  zeros(o+1,1)        zeros(o+1, 1)           -ones(o+1,1)];

        A4 = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)'  -ones(m + 1,1)      zeros(m+1,1)            zeros(m+1,1)];
        A5 = [zeros(n+1,1) -1.*ones(n+1,1) (0:1:n)'  zeros(n + 1,1)      ones(n+1,1)             zeros(n+1,1)];
        A6 = [zeros(o+1,1) -1.*ones(o+1,1) (0:1:o)'  zeros(o + 1,1)      zeros(o+1,1)            ones(o+1,1)];


    case 'lambda_rho'
        f = [1 -1 0  0 0 ];
        A1 = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)'  -ones(m + 1, 1)              zeros(m + 1, 1) ];
        A2 = [ones(n+1,1) zeros(n+1,1) -1.*(0:1:n)'  zeros(n + 1, 1)              zeros(n + 1, 1)];
        A3 = [ones(o+1,1) zeros(o+1,1) -1.*(0:1:o)'  zeros(o+1,1)                   -ones(o+1,1)];

        A4 = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)'  -ones(m + 1,1)                  zeros(m+1,1)];
        A5 = [zeros(n+1,1) -1.*ones(n+1,1) (0:1:n)'  zeros(n + 1,1)                   zeros(n+1,1)];
        A6 = [zeros(o+1,1) -1.*ones(o+1,1) (0:1:o)'  zeros(o + 1,1)                  ones(o+1,1)];
        
        
    case 'lambda_mu'
        f = [1 -1 0  0 0 ];
        A1 = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)'  -ones(m + 1, 1)     zeros(m + 1, 1)        ];
        A2 = [ones(n+1,1) zeros(n+1,1) -1.*(0:1:n)'  zeros(n + 1, 1)     -ones(n + 1, 1)        ];
        A3 = [ones(o+1,1) zeros(o+1,1) -1.*(0:1:o)'  zeros(o+1,1)        zeros(o+1, 1)          ];

        A4 = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)'  -ones(m + 1,1)      zeros(m+1,1)           ];
        A5 = [zeros(n+1,1) -1.*ones(n+1,1) (0:1:n)'  zeros(n + 1,1)      ones(n+1,1)            ];
        A6 = [zeros(o+1,1) -1.*ones(o+1,1) (0:1:o)'  zeros(o + 1,1)      zeros(o+1,1)           ];
        
    case 'mu_rho'
        f = [1 -1 0  0 0 ];
        A1 = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)'       zeros(m + 1, 1)         zeros(m + 1, 1) ];
        A2 = [ones(n+1,1) zeros(n+1,1) -1.*(0:1:n)'       -ones(n + 1, 1)         zeros(n + 1, 1)];
        A3 = [ones(o+1,1) zeros(o+1,1) -1.*(0:1:o)'       zeros(o+1, 1)           -ones(o+1,1)];

        A4 = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)'       zeros(m+1,1)            zeros(m+1,1)];
        A5 = [zeros(n+1,1) -1.*ones(n+1,1) (0:1:n)'       ones(n+1,1)             zeros(n+1,1)];
        A6 = [zeros(o+1,1) -1.*ones(o+1,1) (0:1:o)'       zeros(o+1,1)            ones(o+1,1)];
        
end
        
B1 = log10(abs(fx));
B2 = log10(abs(gx));
B3 = log10(abs(hx));

B4 = log10(abs(fx));
B5 = log10(abs(gx));
B6 = log10(abs(hx));


% % 
% Must first remove any zeros from lambda_vec, and corresponding rows in 
% Part one

% Remove zeros
[B1, A1] = RemoveZeros(B1, A1);
[B2, A2] = RemoveZeros(B2, A2);
[B3, A3] = RemoveZeros(B3, A3);

[B4, A4] = RemoveZeros(B4, A4);
[B5, A5] = RemoveZeros(B5, A5);
[B6, A6] = RemoveZeros(B6, A6);




A = [A1; A2; A3; A4; A5; A6];
b = [ B1; B2; B3; -B4; -B5; -B6];


x = linprog(f, -A, -b);



try
    
    
    switch scaling_method
        
        case 'lambda_mu_rho'
    
    theta = 10^x(3);
    lambda = 10^x(4);
    mu = 10^x(5);
    rho = 10^x(6);
    
        case 'lambda_rho'
            
            theta = 10^x(3);
            lambda = 10^x(4);
            mu = 1;
            rho = 10^x(5);
    
        case 'lambda_mu'
            
            theta = 10^x(3);
            lambda = 10^x(4);
            mu = 10^x(5);
            rho = 1;

        case 'mu_rho'
            
    theta = 10^x(3);
    lambda = 1;
    mu = 10^x(4);
    rho = 10^x(5);
        otherwise
            error('err')
    end
    
catch
    fprintf('Failed to optimize\n')
    theta = 1;
    lambda = 1; 
    mu = 1;
    rho = 1;

end

end


function [f_vec_a, Part1] = RemoveZeros(f_vec_a, Part1)

% Get the index of the zero rows in lambda
index_zero_f_vec_a = find(f_vec_a==0);
% Remove the corresponding zeros from lambda
f_vec_a(index_zero_f_vec_a,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part1(index_zero_f_vec_a,:) = [];

end