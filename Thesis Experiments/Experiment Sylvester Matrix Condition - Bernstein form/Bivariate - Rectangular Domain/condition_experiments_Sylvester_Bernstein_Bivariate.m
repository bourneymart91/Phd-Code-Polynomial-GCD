function [] = condition_experiments_Sylvester_Bernstein_Bivariate(m1, m2, n1, n2)
% Get the condition of a set of Sylvester matrices
%
% %
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% n1 : (Int) Degree of g(x,y) with respect to x
%
% n2 : (Int) Degree of g(x,y) with respect to y
%
% Example
% >> condition_experiments_Sylvester_Bernstein_Bivariate(4,4,3,3)

% Basis : Bernstein Rectangular Domain
% Type : Bivariate


warning('off')

% Set coefficients of f(x,y) and g(x,y) to all be ones
fxy = ones(m1 + 1, m2 + 1);
gxy = ones(n1 + 1, n2 + 1);


arr_Methods = {...
    'DTQ',...
    'DT',...
    'T',...
    'DTQ Denominator Removed'...
    };

% Get number of methods
nMethods = length(arr_Methods);

% Initialise array to store matrices of condition numbers
arr_condition = cell(nMethods);

for i = 1:1:nMethods
    
    method_name = arr_Methods{i};
    
    cond_matrix = zeros(min(m1,n1) , min(m2,n2));
    
    for k1 = 1:1:min(m1,n1)
        for k2 = 1:1:min(m2,n2)

            Sk = BuildSylvesterMatrix(fxy, gxy, k1, k2, method_name);
            cond_matrix(k1,k2) = cond(Sk);
            
        end
        
    end
    
    arr_condition{i} = cond_matrix;
    
end

figure()
hold on
for i = 1:1:nMethods
    
    method_name = arr_Methods{i};
    condition_matrix = arr_condition{i};
    s1 = surf(log10(condition_matrix));
    set(s1,'FaceAlpha',0.9,'DisplayName',method_name);
    
end

legend(gca,'show');
title('Condition of subresultant S_{k_{1},k_{2}}')
xlabel('k1')
ylabel('k2')
zlabel('log_{10} condition')
hold off
grid on
end



function Sk = BuildSylvesterMatrix(fxy, gxy, k1, k2, method_name)

[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);


switch method_name
    case 'DTQ'
        
        % Build the matrix D
        D = BuildD(m1, m2, n1-k1, n2-k2);

        % Build the matrix T
        T1 = BuildT1(fxy, n1-k1, n2-k2);
        T2 = BuildT1(gxy, m1-k1, m2-k2);

        % Build the Matrix Q
        Q = BuildQ(n1-k1,n2-k2,m1-k1,m2-k2);
        
        Sk = D*[T1 T2]*Q;
        
    case 'DT'
        % Build the matrix D
        D = BuildD(m1, m2, n1-k1, n2-k2);

        % Build the matrix T
        T1 = BuildT1(fxy, n1-k1, n2-k2);
        T2 = BuildT1(gxy, m1-k1, m2-k2);

        Sk = D*[T1 T2];
        
        
    case 'T'

        % Build the matrix T
        T1 = BuildT1(fxy, n1-k1, n2-k2);
        T2 = BuildT1(gxy, m1-k1, m2-k2);

        Sk = [T1 T2];
        
    case 'DTQ Denominator Removed'
        % Build the matrix D
        D = BuildD(m1, m2, n1-k1, n2-k2);

        % Build the matrix T
        T1 = BuildT1(fxy, n1-k1, n2-k2);
        T2 = BuildT1(gxy, m1-k1, m2-k2);

        % Build the Matrix Q
        Q = BuildQ(n1-k1,n2-k2,m1-k1,m2-k2);
        
        % With denominators removed
        com_denom_T1 = nchoosek(m1+n1-k1,n1-k1) * nchoosek(m2+n2-k2,n2-k2);
        com_denom_T2 = nchoosek(m1+n1-k1,m1-k1) * nchoosek(m2+n2-k2,m2-k2);
        T1 = T1 ./ com_denom_T1;
        T2 = T2 ./ com_denom_T2;
        Sk = D*[T1 T2]*Q;
        
    otherwise
        error('err')
        
end

end



function D = BuildD(m1,m2,n1_k1,n2_k2)
% Build the matrix D

% Divide each row of d by the corresponding \binom{m1+n1-k1}{i}
% Divide each col of d by the corresponding \binom{m2+n2-k2}{j}
D_mat = GetWithBinomials(ones(m1+n1_k1+1,m2+n2_k2+1));

% Get the D_matrix as a vector
D_vec = GetAsVector(1./D_mat);

% Form a diagonal matrix from the vector D_vec
D = diag(D_vec);

end

function Q1 = BuildQ1(n1_t1,n2_t2)
% Build the matrix Q1, a partition of the matrix Q, where Q is the diagonal
% matrix of binomial coefficients of the vector [f;g].

Q_mat = GetWithBinomials(ones(n1_t1+1,n2_t2+1));

% Get the matrix Q as a vector
Q_vec = GetAsVector(Q_mat);

% Diagonalise Q
Q1 = diag(Q_vec);

end

function T1 = BuildT1(fxy, n1_k1, n2_k2)
% Build the matrix T1 where T1(f) * [v] gives a vector of coefficients of
% the polynomial multiplication f(x,y) * v(x,y).
%
% Inputs.
%
% fww_matrix : Coefficients of polynomial f(w,w)
%
% n1_k1 : Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : Degree of polynomial v(x,y) with repsect to y


% get size of f(x,y)
[m1,m2] = GetDegree_Bivariate(fxy);

% % Build the matrix f(w,w) including binomials.
fxy_bi = GetWithBinomials(fxy);

% Initalise a zero matrix
zero_matrix = zeros(m1 + n1_k1 + 1, m2 + n2_k2 + 1);

% Get the number of rows in T_{k_{1},k_{2}}(f)
nColsT1 = (n1_k1 + 1) * (n2_k2 + 1);
nRowsT1 = (m1+n1_k1 +1 ) * (m2+n2_k2 +1);

T1 = zeros(nRowsT1,nColsT1);

% for every diagonal of the matrix vxy.
nDiags_vxy = (n2_k2+1)+(n1_k1+1) -1;

% Initialise a count
count = 1;

for tot0 = 0:1:nDiags_vxy - 1;
    for i = tot0:-1:0
        j = tot0-i;
        % if j1 is within the number of columns of vxy_mtrx.
        % if j2 is within the number of rows of vxy_mtrx.
        if j <= n2_k2  && i <= n1_k1
            
            % Padd the coefficients
            fxy_matrix_bi_padded_new = zero_matrix;
            fxy_matrix_bi_padded_new((i+1):(m1)+(i+1),(j+1):(m2)+(j+1)) = fxy_bi;
            
            temp_vec = GetAsVector(fxy_matrix_bi_padded_new);
            
            % Assign column
            T1(:,count) = temp_vec;
            
            % Increment Counter
            count = count + 1;
        end
        %
    end
end



end

function Q = BuildQ(n1_k1,n2_k2,m1_k1,m2_k2)
% Build the diagonal matrix Q whose entries contain binomial coefficients
% corresponding to the entries of u(x,y) and v(x,y) in the vector [v;-u]

Q1 = BuildQ1(n1_k1,n2_k2);
Q2 = BuildQ1(m1_k1,m2_k2);

Q = blkdiag(Q1,Q2);

end

function f_vec = GetAsVector(fxy)
% Given the polynomial f(x,y) in the Bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.

% Get the degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Initialise a counter
count = 1;

% Initialise the vector of coefficients of f(x,y)
f_vec = zeros((m1+1)*(m2+1),1);

% Get number of diagonals in fxy_matrix
nDiags_fxy = (m1+1)+(m2+1)-1;

% For each diagonal of f(x,y), read into the vector, starting lower left to
% upper right.
for tot = 0:1:nDiags_fxy
    
    for i = tot:-1:0
        j = tot - i;
        
        if(i > m1 || j > m2)
            % Do nothing
        else
            
            % Assign matrix entry to next available vector entry
            f_vec(count) = fxy(i+1,j+1);
            
            % Increment counter
            count = count + 1;
        end
        
    end
    
end


end

function fxy_bi = GetWithBinomials(fxy)
% Given the matrix of coefficients f(x,y), include the binomial
% coefficients in the matrix.
%
% Inputs
%
%
% fxy : Matrix of coefficients of the polynomial f(x,y)
%
% Outputs.
%
%
% fxy_bi : Matrix of coefficients of the polynomial f(x,y) including
%          binomial coefficients, so fxy_bi is the polynomial coefficients in the
%          modified Bernstein Basis.

% Get the degree of f(x,y)
[m1,m2] = GetDegree_Bivariate(fxy);

% Get binomial coefficients.
bi_m1 = GetBinomials(m1);
bi_m2 = GetBinomials(m2);

% Get f(x,y) with binomials included (in modified Bernstein Basis).
fxy_bi = diag(bi_m1) * fxy * diag(bi_m2);

end

function [bi_m] = GetBinomials(m)
% Calculate the binomial coefficients mC0,...,mCm and return as a
% column vector

% Initialise the column vector of binomials
bi_m = zeros(m,1);

% for each entry of the vector of binomials calculate nchoosek
for i = 0:1:m
    bi_m(i+1) = nchoosek(m,i);
end

end

function [m1,m2] = GetDegree_Bivariate(fxy)
% GetDegree(fxy)
%
% Get the degree strucuture m_{1} and m_{2} of the bivariate polynomial
% f(x,y)
%
% Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
%
% % Outputs
%
% m1 : (Int)
%
% m2 : (Int)

% Get the dimensions of the matrix of coefficients of f(x,y)
[r,c] = size(fxy);

% Get the degree with respect to x.
m1 = r - 1;

% Get the degree with respect to y.
m2 = c - 1;

end





