
function [root_mult_array_fx, root_mult_array_gx , root_mult_array_dx] = Examples_GCD_FromRoots(ex_num)


pattern = 'Custom:m=(\d+).n=(\d+).t=(\d+).low=(-?\d+).high=(-?\d+)';

if ~isempty(regexp(ex_num,pattern,'start'))
    
    str = ex_num;
    
    expression_m = regexp(str,'m=(\d+)','tokens');
    m_str = expression_m{1};
    
    expression_n = regexp(str,'n=(\d+)','tokens');
    n_str = expression_n{1};
    
    expression_t = regexp(str,'t=(\d+)','tokens');
    t_str = expression_t{1};
    
    expression_low = regexp(str, 'low=(-?\d+)', 'tokens');
    low_str = expression_low{1};
    
    expression_high = regexp(str, 'high=(-?\d+)','tokens');
    high_str = expression_high{1};
    
    
    
    m = str2double(m_str);
    n = str2double(n_str);
    t = str2double(t_str);
    intvl_low = str2double(low_str);
    intvl_high = str2double(high_str);
    
    [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,n,t,intvl_low, intvl_high);
    
    
    
    
else
    switch ex_num
        
        case 'Example 7.1'
            % From The computation of the degree of an approximate greatest
            % common divisor of two Bernstein Polynomials - Bourne, Winkler
            % & Yi.
            root_mult_array_fx = ...
                [
                0.10    4
                0.30    2
                0.50    2
                0.70    3
                0.80    2
                2.50    3
                -3.40   3
                ];
            
            root_mult_array_gx = ...
                [
                0.10    3
                0.80    2
                0.85    4
                0.90    4
                1.10    3
                ];
        
        case 'Example 7.2'
            % From The computation of the degree of an approximate greatest
            % common divisor of two Bernstein Polynomials - Bourne, Winkler
            % & Yi.
            
            root_mult_array_fx = ...
                [
                0.1     3
                0.56    4
                0.75    3
                0.82    3
                1.37    3
                -0.27   3
                -1.46   2
                ];
            root_mult_array_gx = ...
                [
                0.1     2
                0.56    4
                0.75    3
                0.99    4
                1.37    3
                2.12    3
                1.2     3
                ];
        
        
            
            
            
        case '1'
            root_mult_array_fx = [
                0.2 2
                0.5 1
                0.7 1
                0.9 1
                1.1 1
                ];
            root_mult_array_gx = [
                0.5 1
                0.1 1
                0.9 1
                0.3 1
                ];
            
            
        case '2'
            root_mult_array_fx = [
                0.2 1
                0.4 1
                0.6 1
                0.8 1
                ];
            
            root_mult_array_gx = [0.2 1
                0.3 1];
            
        case '3'
            root_mult_array_fx = [
                0.56 20
                0.75 3
                0.82 3
                0.37 3];
            
            root_mult_array_gx = [
                0.56    20
                0.75    3
                0.99    4
                0.37    3
                0.12    3
                0.20    3
                ];
            
        case '4'
            
            root_mult_array_fx = [0.1    20
                0.5    2];
            
            root_mult_array_gx = [0.1    20
                0.9    1];
            
        case '5'
            
            root_mult_array_fx = [0.1    2
                0.3     2
                0.5    2];
            
            root_mult_array_gx = [0.1    2];
            
            % From Winkler Paper - Methods for the computation of the degree of an
            % approximate greatest common divisor of two inexact bernstein basis
            % polynomials.
        case '6'
            root_mult_array_fx = [
                0.10    3
                0.56    4
                0.75    3
                0.82    3
                1.37    3
                -0.27   3
                1.46    2
                ];
            
            root_mult_array_gx = [
                0.10    2
                0.56    4
                0.75    3
                0.99    4
                1.37    3
                2.12    3
                1.20    3
                ];
            
        case '7'
            root_mult_array_fx = [
                0.23   4
                0.43   3
                0.57   3
                0.92   3
                1.70   3
                ];
            root_mult_array_gx = [
                0.23   4
                0.30   2
                0.77   5
                0.92   2
                1.20   5
                ];
            
        case '8'
            root_mult_array_fx = [0.1  5
                0.56 4
                0.75 3
                0.82 3
                1.37 3];
            root_mult_array_gx = [0.1  5
                0.56 4
                0.75 3
                0.99 4
                1.37 3
                2.12 3
                1.20 3];
            
            % Example 6.2
        case '9'
            root_mult_array_fx = [
                0.14    3
                0.56    3
                0.89    4
                1.45    4
                2.37    3
                -3.61   5
                ];
            root_mult_array_gx = [
                0.14    4
                0.99    1
                2.37    3
                -0.76   2
                -1.24   2
                -3.61   7
                ];
            
        case '10'
            root_mult_array_fx = [
                0.14    3
                0.56    3
                0.89    4
                0.37    3
                ];
            
            root_mult_array_gx = [
                0.14    3
                0.37    3
                0.76   2
                0.24   2
                ];
            
            
            % Example 6.2
        case '11'
            root_mult_array_fx = [
                0.10    3
                0.56    5
                1.40    4
                1.79    3
                2.69    2
                -2.68   3
                ];
            root_mult_array_gx = [
                0.10    4
                0.56    4
                1.79    3
                2.68    2
                2.69    3
                -1.40   2
                ];
        case '12'
            root_mult_array_fx = [
                0.10    7;
                0.50    12;
                ];
            
            root_mult_array_gx = [
                0.10    6;
                0.50    11;
                0.99    1;
                ];
        case '13'
            root_mult_array_fx = ...
                [
                0.5     4;
                ];
            root_mult_array_gx = ...
                [
                0.5     3;
                ]
            
        case '14'
            
            intvl_low = -1;
            intvl_high = 1;
            
            t = 5;
            m = 10;
            ex_num = 7;
            [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,ex_num,t,intvl_low, intvl_high);
            
        case 'Custom'
            intvl_low = -1;
            intvl_high = 1;
            
            prompt = 'Enter the degree of Polynomial f(x) :';
            m = input(prompt);
            
            prompt = 'Enter the degree of Polynomial g(x) :';
            ex_num = input(prompt);
            
            prompt = 'Enter the degree of Polynomial d(x) :';
            t = input(prompt);
            
            
            %         t = 5;
            %         m = 10;
            %         n = 7;
            [root_mult_array_fx,root_mult_array_gx] = BuildRandomPolynomials(m,ex_num,t,intvl_low, intvl_high)
            
            
        otherwise
            error('error: Example Number not valid')
            
    end
    
end
root_mult_array_dx = GetGCDRoots(root_mult_array_fx,root_mult_array_gx)
root_mult_array_ux = GetQuotientRoots(root_mult_array_fx,root_mult_array_dx);
root_mult_array_vx = GetQuotientRoots(root_mult_array_gx,root_mult_array_dx);

m = sum(root_mult_array_fx(:,2));
ex_num = sum(root_mult_array_gx(:,2));
d = sum(root_mult_array_dx(:,2));
t = d;

end



function d_roots = GetGCDRoots(f_roots,g_roots)
% getDivisor(f_roots,g_roots)
%
% Given the set of roots of f(x) and roots of g(x) and their multiplicities,
% get the common roots of polynomials f(x) and g(x) given by d(x).

% get the number of roots in polynomial f
num_roots_f = size(f_roots,1);

% Initialise the set of roots of d(x)
d_roots = [];

% for each root in f(x), check to see if it exists in g(x)
for i = 1:1:num_roots_f
    
    % Get the root of f(x)
    root = f_roots(i,1);
    
    % Get the multiplicity of root r_{i}
    mult_root_in_f = f_roots(i,2);
    
    % Get the number of distinct roots in g
    [distinct_roots_g,~] = size(g_roots);
    if  distinct_roots_g == 0
        return
    end
    
    % Look if the root r_{i} exists in g(x)
    if ~isempty(find(g_roots(:,1) == root));
        
        % Get the index of the row which corresponds to the root r_{i} in
        % the matrix of roots of g(x)
        [row_d,~] = find(g_roots(:,1) == root);
        
        % Get the multiplicity of the root r_{i} in g(x)
        mult_root_in_g = g_roots(row_d,2);
        
        % Calculate the multiplicity of the root in d(x)
        mult_root_in_d = min(mult_root_in_f,mult_root_in_g);
        
        % Add the root to d(x)
        d_roots = [d_roots ; root mult_root_in_d];
    end
end


end


function u_roots = GetQuotientRoots(f_roots,d_roots)
% Fet the roots of quotient polynomial u(x) given the roots of polynomial
% f(x), and the roots of polynomial d(x), where d(x) is the GCD of f(x) and
% g(x) and where f(x)/u(x) = d(x)

% Get the number of distinct roots in f(x)
nRoots_f_x = size(f_roots,1);

% Initialise an empty matrix of roots in u(x)
u_roots = [];

% Catch the case that the degree of the GCD is zero, and therefore the quotient
% polynomial u(x) is equal to f(x)
[nRoots_d,~] = size(d_roots);
if nRoots_d == 0 % If d has no roots, then d is scalar and u(x) = f(x)
    u_roots = f_roots;
    return
end

% For each root of f(x), check to see if it exists in d(x).
for i = 1:1:nRoots_f_x
    
    % Get the root r_{i}
    root = f_roots(i,1);
    
    % Get the multiplicity of the root r_{i}
    mult_f = f_roots(i,2);
    
    % Look up the root in roots of d
    if ~isempty(find(d_roots(:,1) == root));
        
        % Get the row on which we find the root
        [row_d,~] = find(d_roots(:,1) == root);
        
        % Get the multiplicity of the root
        mult_d = d_roots(row_d,2);
        
        % Subtract multiplicty in d to obtain multiplicity in quotient
        % polynomial u
        mult_u = mult_f - mult_d;
        
        % Add the root and its multiplicity to the set of roots for
        % quotient polynomial u
        if mult_u > 0
            u_roots = [u_roots; root mult_u];
        end
        
    else
        
        u_roots = [u_roots; root mult_f];
    end
end
end
