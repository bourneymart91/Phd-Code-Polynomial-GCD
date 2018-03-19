function [alpha,theta] =  OptimalAlphaTheta(F_max,F_min,G_max,G_min)
% OptimalAlphaTheta(F_max,F_min,G_max,G_min)
%
% This function computes the optimal values alpha and theta for the
% preprocessing operations on the kth Sylvester subresultant matrix
% S_{k}(f,g)Q_{k} = D^{-1}_{k}T_{k}Q_{k}
%
%
% Inputs
%
%
% F_max   :  A vector of length m+1, such that F_max(i) stores the
%            element of maximum magnitude of S(f,g)Q that contains the
%            coefficient a(i) of f, i=1,...,m+1.
%
% F_min   :  A vector of length m+1, such that F_min(i) stores the
%            element of minimum magnitude of S(f,g)Q that contains the
%            coefficient a(i) of f, i=1,...,m+1.
%
% G_max   :  A vector of length n+1, such that G_max(i) stores the
%            element of maximum magnitude of S(f,g)Q that contains the
%            coefficient b(i) of g, i=1,...,n+1.
%
% G_min   :  A vector of length n+1, such that G_min(i) stores the
%            element of minimum magnitude of S(f,g)Q that contains the
%            coefficient b(i) of g, i=1,...,n+1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define vector f
f = [1, -1, 0, 0];

% get degree of polynomials f and g
m = GetDegree(F_max);
n = GetDegree(G_max);

% Assemble the four submatrices.
Part1 = [ones(m+1,1)    zeros(m+1,1)    -(0:1:m)'  zeros(m+1,1)];
Part2 = [ones(n+1,1)    zeros(n+1,1)    -(0:1:n)'  -ones(n+1,1)];
Part3 = [zeros(m+1,1)   -ones(m+1,1)     (0:1:m)'  zeros(m+1,1)];
Part4 = [zeros(n+1,1)   -ones(n+1,1)     (0:1:n)'   ones(n+1,1)];

A = -[Part1;Part2;Part3;Part4];

% Construct the vector b.
b= -[
     log10(abs(F_max));
     log10(abs(G_max));
    -log10(abs(F_min));
    -log10(abs(G_min))];

% Solve the linear programming problem and extract alpha and theta
% from the solution vector x.

    warning('off')


    x = linprog(f,A,b);
    warning('on')
    
    if numel(x) == 4
        
        theta = real(10^x(3));
        alpha = real(10^x(4));
        
    elseif numel(x) <4
        % Exclude alpha, optimise theta
        Part1 = [ ones(m + 1, 1)  zeros(m + 1, 1)  (0 : -1 : -m)'   ];
        Part2 = [ ones(n + 1, 1)  zeros(n + 1, 1)  (0 : -1 : -n)'   ];
        Part3 = [ zeros(m + 1, 1) -ones(m + 1, 1)  (0 : 1 : m)'     ];
        Part4 = [ zeros(n + 1, 1) -ones(n + 1, 1)  (0 : 1 : n)'     ];
        

        
        f = [1 -1 0];
        A = -[Part1; Part2; Part3; Part4];
        x = linprog(f,A,b);
        
        if numel(x) == 3
            theta = real(10^x(3));
            alpha = 1;
            
        elseif numel(x) <3
            % Exclude theta, optimise alpha
            Part1 = [ ones(m+1,1)   zeros(m+1,1) zeros(m+1,1)];
            Part2 = [ ones(n+1,1)   zeros(n+1,1) -ones(n+1,1)];
            Part3 = [ zeros(m+1,1)  ones(m+1,1)  zeros(m+1,1)];
            Part4 = [ zeros(n+1,1)  ones(n+1,1)  ones(n+1,1)];
            
            f = [1 -1 0];
            A = -[Part1; Part2; Part3; Part4];
            x = linprog(f,A,b);
            
            if numel(x) == 3
                alpha = real(10^x(3));
                theta = 1;
            elseif numel(x) <3
                alpha = 1;
                theta = 1;
            end
            
            
        end
        
    end




end