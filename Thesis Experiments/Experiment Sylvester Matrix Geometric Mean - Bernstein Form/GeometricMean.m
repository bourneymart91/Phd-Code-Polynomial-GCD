function [] = GeometricMean(m,n)

vMu = zeros(min(m,n),1);
vMu_matlab = zeros(min(m,n),1);

vGM_coefficients = zeros(m+1,min(m,n));




for k = 1:1:min(m,n)

    vMu(k) = GetGeometricMeanofCf_log(m,n-k);
    vMu_matlab(k) = GetGeometricMeanOfCf_matlab(m,n-k);
    
    % Return the vector of the geometric mean of each of the coefficients
    % a_{i}.
    vGM_coefficients(:,k) = GetGeometricMean_Coefficients(m,n-k);
    
    
end


figure()
hold on
plot(log10(vMu_matlab),'-s')
hold off


figure()
hold on
plot(log10(vGM_coefficients'))
hold off

end

function mu = GetGeometricMeanOfCf_matlab(m,n_k)

Matrix = zeros(m+1,n_k+1);

% for each coefficient a_{0},...,a_{m}
for i  = 0 :1: m
    % for each column 0,...,n-k
    for j = 0:1:n_k
        Matrix(i+1,j+1) = (nchoosek(m,i) ...
            * nchoosek(n_k,j)) ...
            / nchoosek(m+n_k,i+j);
    end
end

% for each coefficient a_{0},...,a_{m}
for i  = 0 :1: m
    % for each column 0,...,n-k
    for j = 0:1:n_k
        Matrix2(i+1,j+1) = (nchoosek(i+j,i) ...
            * nchoosek(m+n_k-(i+j),m-i)) ...
            / nchoosek(m+n_k,m);
    end
end


vec = reshape(Matrix,1,(m+1)*(n_k+1));

mu = geomean(vec);


end

function [mu] = GetGeometricMeanofCf_log(m,n_k)

part_a = 0;
for i = 0:1:m
    for j = 0:1:n_k
        part_a = part_a + lnnchoosek(i+j,j);
    end
end

part_a_log = (1./((n_k+1)*(m+1)))* (part_a);

part_b = 0;
for i = 0:1:m
    for j = 0:1:n_k
        part_b = part_b + lnnchoosek(m+n_k-(i+j),m-i);
    end
end
part_b_log = (1./((n_k+1)*(m+1)))* (part_b);


part_c_log = lnnchoosek(m+n_k,m);


ln_mu = part_a_log + part_b_log - part_c_log;


mu = 10.^(ln_mu);


end

function nCk = lnnchoosek(n,k)
nCk = sum(log10(1:1:n)) - sum(log10(1:1:k)) - sum(log10(1:1:(n-k)));
end


function [vGM] = GetGeometricMean_Coefficients(m,n_k)
% Get the geometric mean of each of the coefficients.

% for each coefficient get the geometric mean of its entries in c(f)
% for each column 0,...,n-k
% for each coefficient a_{0},...,a_{m}
for i  = 0 :1: m
    % for each column 0,...,n-k
    for j = 0:1:n_k
        Matrix(i+1,j+1) = (nchoosek(m,i) ...
            * nchoosek(n_k,j)) ...
            / nchoosek(m+n_k,i+j);
    end
end

vGM = geomean(Matrix,2);

end


