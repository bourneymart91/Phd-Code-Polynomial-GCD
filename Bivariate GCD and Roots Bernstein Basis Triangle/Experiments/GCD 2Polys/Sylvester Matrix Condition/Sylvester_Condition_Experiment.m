function [] = Sylvester_Condition_Experiment()
% Sylvester_Condition_Experiment(m,n)
%
% % Inputs
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% Examples
%
% >> Sylvester_Condition_Experiment(4,3)

% Basis : Bernstein
% Type : Bivariate Triangle

warning('off')

% Add examples folder and all subfolders in current folder
addpath(genpath('../Examples'))
addpath(genpath(pwd));

% % % Get Example
% ex_num = '1';
% [f_root_mult_arr,g_root_mult_arr,d_root_mult_arr,...
%     u_root_mult_arr,v_root_mult_arr] = GCD_Examples_Bivariate_2Polys(ex_num);
% 
% % Get coefficients of the polynomials
% [fxy] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr);
% [gxy] = GetCoefficientsFromSymbolicRoots(g_root_mult_arr);
% [dxy] = GetCoefficientsFromSymbolicRoots(d_root_mult_arr);
% [uxy] = GetCoefficientsFromSymbolicRoots(u_root_mult_arr);
% [vxy] = GetCoefficientsFromSymbolicRoots(v_root_mult_arr);
% 
% % Get the symbolic polynomials in power form
% symbolic_f = GetSymbolicPoly(f_root_mult_arr);
% symbolic_g = GetSymbolicPoly(g_root_mult_arr);
% symbolic_d = GetSymbolicPoly(d_root_mult_arr);
% display(symbolic_f)
% display(symbolic_g)
% display(symbolic_d)
% 
% 
% m = GetDegree_Bivariate(fxy);
% n = GetDegree_Bivariate(gxy);
% 
% 
% v_fxy = GetAsVector(fxy);
% v_gxy = GetAsVector(gxy);

m = 15;
n = 30;


% Set coefficients of f(x,y) and g(x,y) to all be ones
nCoefficients_fxy = nchoosek(m + 2, 2);
nCoefficients_gxy = nchoosek(n + 2, 2);

%v_fxy = 10.*rand(nCoefficients_fxy, 1);
%v_gxy = 2.*rand(nCoefficients_gxy, 1);

v_fxy = ones(nCoefficients_fxy, 1);
v_gxy = ones(nCoefficients_gxy, 1);


% Put coefficients of f(x,y) into a matrix of dimension m+1 x m+1
temp_vec = ...
    [
        v_fxy;
        zeros(nchoosek(m+1,2),1);
    ];
mat_fxy = GetAsMatrix(temp_vec,m,m);

% Put coefficients of g(x,y) into a matrix of dimension n+1 x n+1
temp_vec = ...
    [
        v_gxy;
        zeros(nchoosek(n+1,2),1);
    ];
mat_gxy = GetAsMatrix(temp_vec,n,n);

% %
% 


arrSylvesterFormat = {'DTQ', 'DT', 'TQ', 'T', 'DTQ Denominator Removed'};
nSylvesterFormats = length(arrSylvesterFormat);

arrConditionSk = cell(nSylvesterFormats, 1);
arrConditionCk = cell(nSylvesterFormats, 1);




for i = 1:1:nSylvesterFormats

    for k = 1:1:min(m,n)

        format_name = arrSylvesterFormat{i};
        arrConditionSk{i}(k) = GetConditionSylvesterMatrix(mat_fxy, mat_gxy, m, n, k, format_name);
        arrConditionCk{i}(k) = GetConditionFirstPartition(mat_fxy,m, n,k, format_name);
        

    end

end

% Plot condition of the Sylvester matrices
figure('name', [mfilename ' : ' 'Conditions'])
hold on
for i = 1:1:nSylvesterFormats
    
    
    vConditionSk = arrConditionSk{i};
    format_name = arrSylvesterFormat{i};
    
    plot(log10(vConditionSk), '-s', 'DisplayName', format_name);
    
end
legend(gca,'show');

hold off




% Plot Conditon of the first partition
figure('name', [mfilename ' : ' 'Conditions'])
hold on
for i = 1:1:nSylvesterFormats
    
    
    vConditionCk = arrConditionCk{i};
    format_name = arrSylvesterFormat{i};
    
    plot(log10(vConditionCk), '-s', 'DisplayName', format_name);
end
legend(gca,'show');

hold off

end

function condition = GetConditionFirstPartition(mat_fxy, m, n, k, sylvester_matrix_variant)
%
% % Inputs
%
% mat_fxy : (Matrix) 
%
% n_k : (Int)
%
% sylvester_matrix_variant: (String)

switch sylvester_matrix_variant

    
    case 'DTQ'
        
        D = BuildD(m, n-k);

        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        
        % Build the matrix Q
        Q1 = BuildQ1(n - k);

        Cf = D*[T1]*Q1;
        
        condition = cond(Cf);
        
    case 'DT'

        D = BuildD(m, n-k);

        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
   
        Cf = D*[T1];
        
        condition = cond(Cf);
        
    case 'TQ'
        
        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        
        % Build the matrix Q
        Q = BuildQ1(n-k);

        Cf = [T1]*Q;
        
        condition = cond(Cf);
        
    case 'T'
        
        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        Cf = T1 ;
        
        condition = cond(Cf);
        
    case 'DTQ Denominator Removed'
        
        % Build the matrix D^{-1}_{m+n-k}
        D = BuildD(m, n - k);

        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);

        % Build the matrix Q
        Q1 = BuildQ1(n-k);
        
        
        % With denominators removed
        com_denom_T1 = nchoosek(m + n - k , n - k);
        
        T1 = T1 ./ com_denom_T1;
        Cf = D*[T1]*Q1;
        condition = cond(Cf);
    otherwise
        error('err')
end

end


function condition = GetConditionSylvesterMatrix(mat_fxy, mat_gxy, m ,n , k, sylvester_matrix_variant)

switch sylvester_matrix_variant

    
    case 'DTQ'
        
        D = BuildD(m, n-k);

        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        T2 = BuildT1(mat_gxy, n, m - k);

        % Build the matrix Q
        Q = BuildQ(m, n, k);

        Sk = D*[T1 T2]*Q;
        
        condition = cond(Sk);
        
    case 'DT'

        D = BuildD(m, n-k);

        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        T2 = BuildT1(mat_gxy, n, m - k);

        Sk = D*[T1 T2];
        
        condition = cond(Sk);
        
    case 'TQ'
        
        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        T2 = BuildT1(mat_gxy, n, m - k);

        % Build the matrix Q
        Q = BuildQ(m, n, k);

        Sk = [T1 T2]*Q;
        
        condition = cond(Sk);
        
    case 'T'
        
        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        T2 = BuildT1(mat_gxy, n, m - k);

        Sk = [T1 T2];
        
        condition = cond(Sk);
        
    case 'DTQ Denominator Removed'
        
        % Build the matrix D^{-1}_{m+n-k}
        D = BuildD(m, n-k);

        % Build the matrix T
        T1 = BuildT1(mat_fxy, m, n - k);
        T2 = BuildT1(mat_gxy, n, m - k);

        % Build the matrix Q
        Q = BuildQ(m, n, k);
        
        
        % With denominators removed
        com_denom_T1 = nchoosek(m + n - k , n - k);
        com_denom_T2 = nchoosek(m + n - k, m - k);
        T1 = T1 ./ com_denom_T1;
        T2 = T2 ./ com_denom_T2;
        Sk = D*[T1 T2]*Q;
        condition = cond(Sk);
        
    otherwise 
        error('err');
end

end



function D = BuildD(m, n)
% Build the diagonal matrix D^{-1} for the convolution of two polynomials
% f(x,y) and g(x,y) of degrees m and n.

temp_mat = zeros(m + n + 1, m + n + 1);

for i = 0 : 1 : m + n
    
    for j = 0 : 1 : m + n - i
    
        temp_mat(i+1,j+1) = 1./Trinomial(m+n,i,j);
    end
    
end


vTrinomials = GetAsVector(temp_mat);

% Only want the nchoosek(m+n+2,2) non-zero terms
nNonZeroTerms = nchoosek(m + n + 2, 2);

vTrinomials = vTrinomials(1 : nNonZeroTerms);

D = diag(vTrinomials);

end

function Q = BuildQ(m, n, k)
%
% % Inputs
%
% m : (Int) Degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% k : (Int) 

Q1 = BuildQ1(n - k);

Q2 = BuildQ1(m - k);

% Build matrix \hat{Q}_{k}
Q = blkdiag(Q1,Q2);

end

function Q = BuildQ1(m)
% Given the coefficients of polynomial f(x,y), build the diagonal matrix Q
% of trinomial coefficients.
%
% % Inputs
%
% m : (Int) Degree of f(x,y)

temp_mat = zeros(m + 1, m + 1);

for i = 0 : 1 : m
    for j = 0 : 1 : m - i
    
        temp_mat(i + 1, j + 1) = Trinomial(m, i, j);
    
    end
end

vect = GetAsVector(temp_mat);

% Get number of nonzeros
nNonZeros = nchoosek(m + 2, 2);

%
vect = vect(1 : nNonZeros);

Q = diag(vect);


end

function Cf = BuildT1(fxy, m, n)
% Build the matrix C where C(f(x,y))*g = h. 

% Inputs.
%
% fxy : (Matrix) coefficients of polynomial f(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)

% Get number of coefficients of nchoosek
nCoefficients_gxy = nchoosek(n + 2, 2);

% Get number of coefficients in the product fg = h(x,y)
nCoefficients_hxy = nchoosek(m + n + 2, 2);

% Initialise a zero matrix
zero_mat = zeros(m + n + 1, m + n + 1);

Cf = zeros(nCoefficients_hxy, nCoefficients_gxy);

% Get fxy with trinomial coefficients
fxy_tri = GetWithTrinomials(fxy,m);


count = 1;

% For each diagonal of coefficients
for k = 0 : 1 : n 
    
    % i : Row index
    % j : Col index
    for i = k : -1 : 0
        j = k - i;
        
        % Get matrix of coefficients of f(x,y) shifted down by i rows and
        % across by j columns
        temp_mat = zero_mat;
        temp_mat(i + 1 : i + m + 1, j + 1 : j + m + 1) = fxy_tri;
        
        temp_vec = GetAsVector(temp_mat);
        
        % Remove all but the first nchoosek(m+n+2,2) coefficients
        temp_vec = temp_vec(1 : nCoefficients_hxy);
        
        % Insert coefficients into the i,j th column of C(f(x,y)).
        
        Cf(:, count) = temp_vec;
        
        % Increment counter
        count = count + 1;
        
    end
end



end


function fxy_matrix = GetAsMatrix(f_vec, m1, m2)
% Given the vector of coefficients of the polynomial f(x,y), 
% format the coefficients as a matrix.
%
% % Inputs
%
% f_vec : (Vector) Vector of coefficients of f(x,y)
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y

% Initialise an empty matrix fxy
fxy_matrix = zeros(m1 + 1, m2 + 1);

% Intialise a counter which will go through each entry of f_vec (The vector
% of coefficients of of f).
count = 1;

% get number of diagonals in the matrix fxy.
nDiagonals_fxy = (m1 + 1) + (m2 + 1) -1;

for k = 0:1:nDiagonals_fxy - 1
    
    for i = k : -1 : 0
        j = k - i;
        
        if i > m1 || j > m2
            % restrict to only the i and j values within the matrix.
        else
            
            fxy_matrix(i + 1, j + 1) = f_vec(count);
            count = count + 1;
        end
    end
end



end

function f_vec = GetAsVector(fxy_matrix)
% Given the polynomial f in the bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.
%
% % Inputs
%
% fxy_matrix : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs
%
% f_vec : (Vector) Coefficients of f(x,y) in vector form


% Get degree of polynomial f(x,y) with respect to x
[m1, m2] = GetDegree_Bivariate(fxy_matrix);

% Initialise a count
count = 1;

% Initialise an empty vector for coefficients of f(x,y)
f_vec = zeros((m1+1)*(m2+1),1);

% Get the number of diagonals in the matrix of coefficients of f(x,y)
nDiagonals = (m1+1)+(m2+1)-1;


for tot = 0 : 1 : nDiagonals
    for i = tot : -1 : 0
        j = tot - i;
        
        if(i > m1 || j > m2)
        else
            f_vec(count) = fxy_matrix(i+1,j+1);
            count = count + 1;
        end
        
    end
    
end


end

function nCk = Trinomial(m, i, j)

nCk = factorial(m)./(factorial(i)*factorial(j)*factorial(m-i-j));

end




function fxy_tri = GetWithTrinomials(fxy,m)
% Given the coefficients of the polynomial f(x,y) in Bernstein form, get
% the coefficients in the scaled Bernstein form. That is, coefficients of
% f(x,y) with trinomials included.

for i = 0 : 1 : m
    for j = 0 : 1 : m - i
        
        fxy_tri(i + 1, j + 1) = fxy(i + 1, j + 1) .* Trinomial(m, i, j);
        
    end
end

end
