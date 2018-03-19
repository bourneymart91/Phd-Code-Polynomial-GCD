function [] = condition_experiments(m,n)
warning('off')

% Preassign vectors
cond_vecD = zeros(1,min(m,n));
cond_vecQ = zeros(1,min(m,n));
cond_vec_S_old_WithDenom = zeros(1,min(m,n));
cond_vec_S_new_WithoutDenom = zeros(1,min(m,n));
cond_vec_S_new_WithDenom = zeros(1,min(m,n));
cond_vec_S_withoutQ = zeros(1,min(m,n));

% Get the condition of the matrix D^{-1} only. For each value
% k=1,...,min(m,n).
for k = 1:1:min(m,n)
   cond_vecD(k) = condition_experimentD(m,n,k); 
end

% Get the condition of the matrix Q only. For each value
% k=1,...,min(m,n). Where the matrix Q is a diagonal matrix and consists of
% the binomial coefficients corresponding to u_{k} and v_{k}.
for k = 1:1:min(m,n)
    cond_vecQ(k) = condition_experimentQ(m,n,k);
end

% Get the condition of the Sylvester subresultant matrices S_{k}(f,g) for k
% = 1,...,min(m,n). Where the entries of the subresultant matrices have
% binomial coefficients in the original format a_{i-j} (m,i-j) (n-k,j) /
% (m+n-k,i)
% S = D^{-1}T(f,g)Q_{k}
for k = 1:1:min(m,n)
    SylvesterMatrix_Old_WithDenom = Build_SylvesterMatrix_Old_WithDenom(m,n,k);
    cond_vec_S_old_WithDenom(k) = cond(SylvesterMatrix_Old_WithDenom);
end

% Get the condition of the Sylvester subresultant matrices S_{k}(f,g) for k
% = 1,...,min(m,n). Where the entries of the subresultant matrices have
% binomial coefficients in the new format, such that the entries of each 
% partition have common denominators.
for k = 1:1:min(m,n)
    SylvesterMatrix_new_WithDenom = Build_SylvesterMatrix_New_WithDenom(m,n,k);
    cond_vec_S_new_WithDenom(k) = cond(SylvesterMatrix_new_WithDenom) ;
end


% Get the condition of the Sylvester subresultant matrices S_{k}(f,g) for k
% = 1,...,min(m,n). Where the entries of the subresultant matrices have
% binomial coefficients in the new format, such that the entries of each 
% partition have common denominators which have been removed.
for k = 1:1:min(m,n)
    SylvesterMatrix_new_WithoutDenom = Build_SylvesterMatrix_New_WithoutDenom(m,n,k);
    cond_vec_S_new_WithoutDenom(k) = cond(SylvesterMatrix_new_WithoutDenom) ;
end

% Get the condition of the Sylvester subresultant matrices S_{k}(f,g) for k
% = 1,...,min(m,n). Where the entries of each subresultant are in the
% original format, and excludes the binomial coefficients corresponding to
% u_{k} and v_{k} S = D^{-1}T(f,g)
for k = 1:1:min(m,n)
    cond_vec_S_withoutQ(k) = Build_SylvesterMatrix_Old_WithoutQ(m,n,k);
end


% Plot 
figure(1)
hold on
title('Comparing the condition of the Sylvester subresultants formed by new rearrangement, vs the old format')
plot(log10(cond_vec_S_old_WithDenom),'red-s', 'DisplayName','Sylvester Old Format')
plot(log10(cond_vec_S_new_WithDenom),'blue-o','DisplayName','Sylvester New Format with denominator')
plot(log10(cond_vec_S_new_WithoutDenom),'green-s','DisplayName','Sylvester New Format Without Denominator')
xlabel('k')
ylabel('log_{10} Condition')
legend(gca,'show')
hold off

% Plot
figure(2)
hold on
title('Comparing the condition of the Sylvester subresultants in the new format with and without the denominator')
plot(log10(cond_vec_S_new_WithDenom),'red-s', 'DisplayName','Sylvester New Format With Denominator')
plot(log10(cond_vec_S_new_WithoutDenom),'blue-s','DisplayName','Sylvester New Format Without Denominator')
xlabel('k')
ylabel('log_{10} Condition')
legend(gca,'show')
hold off

% Plot
figure(3)
hold on
title('Comparing condition with/without inclusion of Q in Sylvester Matrix')
plot(log10(cond_vec_S_old_WithDenom),'red-s', 'DisplayName','With Q')
plot(log10(cond_vec_S_withoutQ),'blue-s','DisplayName','Without Q')
legend(gca,'show')
xlabel('k')
ylabel('log_{10} Condition')
hold off

end


function [condition] = condition_experimentD(m,n,k)
% Get the condition of matrix D^-1

D = zeros(1,m+n-k+1);

for i = 0:1:m+n-k+1
   D(i+1) = nchoosek(m+n-k+1,i);
end

D = 1./D;

D = diag(D);
condition = cond(D);


end

function [condition] = condition_experimentQ(m,n,k)

Q1 = zeros(1,n-k+1);
Q2 = zeros(1,m-k+1);

for i = 0:1:n-k
    Q1(i+1) = nchoosek(n-k,i) ;
end

for i = 0:1:m-k
    Q2(i+1) = nchoosek(m-k,i);
end

Q = [Q1 Q2];
Q = diag(Q);

condition = cond(Q);
end

function [s] = Build_SylvesterMatrix_Old_WithDenom(m,n,k)
% for each column of first partition

s = zeros(m+n-k+1,m+n-2*k+2);

for j = 0:1:n-k
    for i = j:1:m+j
        a(i+1,j+1) = (nchoosek(m,i-j) * nchoosek(n-k,j)) ./ nchoosek(m+n-k,i);
    end
end

for j = 0:1:m-k
    for i = j:1:n+j
        b(i+1,j+1) = (nchoosek(n,i-j) * nchoosek(m-k,j)) ./ nchoosek(m+n-k,i);
    end
end

s = [a,b];


end


function [s] = Build_SylvesterMatrix_New_WithDenom(m,n,k)

s = zeros(m+n-k+1,m+n-2*k+2);

% for each column of first partition
for j = 0:1:n-k
    for i = j:1:m+j
        a(i+1,j+1) = nchoosek(m+n-k-i,n-k-j) * nchoosek(i,j);
    end
end

for j = 0:1:m-k
    for i = j:1:n+j
        b(i+1,j+1) = nchoosek(m+n-k-i,m-k-j) * nchoosek(i,j);
    end
end

s = [a./nchoosek(m+n-k,n-k),b./nchoosek(m+n-k,m-k)];
%s = [a,b];

end

function [s] = Build_SylvesterMatrix_New_WithoutDenom(m,n,k)
% Build the Sylvester subresultant matrix S_k(f,g), where the coefficients
% of f and g are set to one. This is in the new format, where the binomial
% coefficients have been rearranged, and excluded the denomiators which are
% common to the two partitions.


s = zeros(m+n-k+1,m+n-2*k+2);

% for each column of first partition
for j = 0:1:n-k
    for i = j:1:m+j
        a(i+1,j+1) = nchoosek(m+n-k-i,n-k-j) * nchoosek(i,j);
    end
end

for j = 0:1:m-k
    for i = j:1:n+j
        b(i+1,j+1) = nchoosek(m+n-k-i,m-k-j) * nchoosek(i,j);
    end
end

s = [a,b];


end

function [condition] = Build_SylvesterMatrix_Old_WithoutQ(m,n,k)

s = zeros(m+n-k+1,m+n-2*k+2);

% for each column of first partition
for j = 0:1:n-k
    for i = j:1:m+j
        a(i+1,j+1) = nchoosek(m,i-j) ./ nchoosek(m+n-k,i);
    end
end

for j = 0:1:m-k
    for i = j:1:n+j
        b(i+1,j+1) = nchoosek(n,i-j) ./ nchoosek(m+n-k,i);
    end
end

s = [a,b];

condition = cond(s);
end

