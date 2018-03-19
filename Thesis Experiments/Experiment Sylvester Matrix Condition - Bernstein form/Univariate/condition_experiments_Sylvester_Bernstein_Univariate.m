function [] = condition_experiments_Sylvester_Bernstein_Univariate(ex_num)
% These experiments get the condition of a variety of forms of the
% Sylvester matrix for polynomials f(x) and g(x) in Bernstein form.
%
% % Inputs
%
% m : (Int)  Degree of the polynomial f(x)
%
% n : (Int) Degree of the polynomial g(x)
%
% Example
%
% condition_experiments_Sylvester_Bernstein_Univariate(20,5)

% Basis - Bernstein
% Type - Univariate



addpath(genpath('../../Examples'))

% Get roots and multiplicities from example file
[f_root_mult_arr, g_root_mult_arr, d_root_mult_arr, ~, ~] = ...
    GCD_Examples_Univariate_2Polys(ex_num);

% Get the coefficients of the polynomials f(x), g(x) and d(x).
fx = BuildPolyFromRootsSymbolic(f_root_mult_arr);
gx = BuildPolyFromRootsSymbolic(g_root_mult_arr);

% Get the symbolic polynomials
fx_sym = GetSymbolicPolyFromSymbolicRoots(f_root_mult_arr);
gx_sym = GetSymbolicPolyFromSymbolicRoots(g_root_mult_arr);
dx_sym = GetSymbolicPolyFromSymbolicRoots(d_root_mult_arr);

display(fx_sym)
display(gx_sym)
display(dx_sym)
warning('off')


m = GetDegree(fx);
n = GetDegree(gx);

%fx = ones(m + 1, 1);
%gx = ones(n + 1, 1);


% Get the condition of the Sylvester subresultant matrices S_{k}(f,g) for k
% = 1,...,min(m,n). Where the entries of the subresultant matrices have
% binomial coefficients in the new format, such that the entries of each
% partition have common denominators which have been removed.




arr_Sylvester_Build_Method = ...
    {...
    'DT',...
    'DTQ',...
    'TQ',...
    'T'...
    'DTQ Denominator Removed',...
    'DTQ Rearranged'
    };

nMethods = length(arr_Sylvester_Build_Method);
arr_cond_vec = cell(nMethods,1);

for i = 1:1:nMethods
    
    method_name = arr_Sylvester_Build_Method{i};
    arr_cond_vec{i} = GetCondition(fx, gx, method_name);
    
end


% Plot
figure_name = sprintf(mfilename);
figure('name', figure_name)
hold on
title('Comparing the condition of the Sylvester subresultants formed by new rearrangement, vs the old format')

for i = 1:1:nMethods
    method_name = arr_Sylvester_Build_Method{i};
    plot(log10(arr_cond_vec{i}),'-s','DisplayName',method_name);
end
xlabel('k')
ylabel('log_{10} Condition')
legend(gca,'show')
hold off




end



function [D] = BuildD(m,n_k)
    D = diag(1./GetBinomials(m+n_k));
end

function T1 = BuildT1(fx,n_k)
%
% % Inputs
%
% fx : (Vector) Vector of coefficients of the polynomial f(x)
%
% n_k : (Int) 

% Get degree of f(x)
m = GetDegree(fx);

% Get f(x) with binomial coefficients
f_bi = fx .* GetBinomials(m);


% Build the matrix T_{n-k}(f(x))
T1 = zeros(m+n_k+1,n_k+1);
for j = 1:1:n_k+1
    T1(j:j+m,j) = f_bi;
end

end

function T = BuildT(fx, gx, k)
%
% % Inputs
%
% fx : (Vector)
%
% gx : (Vector) 
%
% k : (Int) 

% Get degree of f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Build the matrix T_{k}(f(x),g(x))
T1 = BuildT1(fx, n - k);
T2 = BuildT1(gx, m - k);
T = [T1 T2];

end

function Q1 = BuildQ1(n_k)
% Build the matrix Q_{1}
%
% % Inputs
%
% n_k : (Int)

Q1 = diag(GetBinomials(n_k));

end

function [binoms] = GetBinomials(m)
% Get the binomial coefficients 
%
% % Inputs
%
% m : (Int)

binoms = zeros(m+1,1);
for i = 0:1:m
    binoms(i+1) = nchoosek(m,i);
end
end




function DT = BuildDT(fx, gx, k)
%
% % Inputs
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% k : (Int)

% Get degree of f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Build the matrix D_{m+n-k}T_{k}(f(x),g(x))
D = BuildD(m, n-k);
T = BuildT(fx, gx, k);
DT = D*T;


end

function TQ = BuildTQ(fx, gx, k)
%
% % Input
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% k : (Int)

% Get the degree of f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Build the matrix T_{k}(f(x),g(x)) \hat{Q}_{k}
T1 = BuildT1(fx,n-k);
T2 = BuildT1(gx,m-k);

Q = BuildQ(m,n,k);

TQ = [T1 T2] * Q;

end


function DTQ = BuildDTQ(fx, gx, k)
%
% % Inputs
%
% fx : (Vector)
%
% gx : (Vector)
%
% k : (Int)

% Get degree of the polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Build the matrix D_{m+n-k}T(f(x),g(x))Q_{k}
D = BuildD(m, n - k);
T1 = BuildT1(fx, n - k);
T2 = BuildT1(gx, m - k);
Q = BuildQ(m, n, k);

DTQ = D*[T1 T2]*Q;

end

function Q = BuildQ(m, n, k)
%
% % Inputs
%
% m : (Int)
%
% n : (Int)
%
% k : (Int)

Q1 = BuildQ1(n - k);
Q2 = BuildQ1(m - k);
Q = blkdiag(Q1, Q2);

end



function [DTQ] = BuildDTQ_new(fx, gx, k)
%
% % Inputs
%
% fx : (Vector)
%
% gx : (Vector)
%
% k : (Int)

m = GetDegree(fx);
n = GetDegree(gx);

a = BuildDT1Q1_new_withDenom(fx, n - k);
b = BuildDT1Q1_new_withDenom(gx, m - k);

DTQ = [a b];
end

function a = BuildDT1Q1_new_withDenom(fx, n_k)
%
% % Inputs
%
% m : (Int)
%
% n_k : (Int)

m = GetDegree(fx);

% for each column of first partition
for j = 0 : 1 : n_k
    for i = j : 1 : m+j
        a(i + 1,j + 1) = fx(i-j+1) * nchoosek(m + n_k - i, n_k - j) * nchoosek(i, j);
    end
end

a = a./ nchoosek(m + n_k, n_k);

end

function [DTQ] = BuildDTQ_new_excludeDenom(fx, gx, k)
% Build the Sylvester subresultant matrix S_k(f,g), where the coefficients
% of f and g are set to one. This is in the new format, where the binomial
% coefficients have been rearranged, and excluded the denomiators which are
% common to the two partitions.

m = GetDegree(fx);
n = GetDegree(gx);

DT1Q1 = BuildDT1Q1_new_excludeDenom(fx, n - k);
DT2Q2 = BuildDT1Q1_new_excludeDenom(gx, m - k);

DTQ = [DT1Q1 , DT2Q2];


end

function DT1Q1 = BuildDT1Q1_new_excludeDenom(fx, n_k)
%
% % Inputs
%
% m : (Int)
%
% n_k : (Int)

m = GetDegree(fx);

% Initialise empty matrix;
DT1Q1 = zeros(m + n_k + 1, n_k + 1);

% for each column of first partition
for j = 0 : 1 : n_k
    for i = j:1:m+j
        DT1Q1(i + 1,j + 1) = fx(i-j+1) .* nchoosek(m + n_k - i, n_k - j) * nchoosek(i,j);
    end
end



end


function m = GetDegree(fx)
%
% % Inputs
% 
% fx : (Int)

m = size(fx,1) - 1;

end


function [vCondition] = GetCondition(fx, gx, method_name)
%
% % Inputs
%
% method_name : (String)
%
% fx : (Vector)
%
% gx : (Vector)
%
% % Outputs
%
% vCondition (Vector) Condition of each S_{k}


m = GetDegree(fx);
n = GetDegree(gx);

vCondition = zeros(min(m,n),1);

switch method_name
    
    case 'DT'
        for k = 1:1:min(m,n)
            
            vCondition(k) = cond(BuildDT(fx, gx, k));
        end
        
    case 'DTQ'
        
        for k = 1:1:min(m,n)
            vCondition(k) = cond(BuildDTQ(fx, gx, k));
        end
        
    case 'TQ'
        for k = 1:1:min(m,n)
            vCondition(k) = cond(BuildTQ(fx, gx, k));
            
        end
        
    case 'T'
        for k = 1:1:min(m,n)
            vCondition(k) = cond(BuildT(fx, gx, k));
            
        end
        
    case 'DTQ Denominator Removed'
        for k = 1:1:min(m,n)
            vCondition(k) = cond(BuildDTQ_new_excludeDenom(fx, gx, k));
            
        end
        
    case 'DTQ Rearranged'
        for k = 1:1:min(m,n)
            vCondition(k) = cond(BuildDTQ_new(fx,gx,k)) ;
        end
    otherwise
        error('err')
end
end

