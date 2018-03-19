function [Y] = BuildDY_Roots(m,n,t,q,x,alpha,theta,ratio)


if q <= n-t+1
    %fprintf('BuildY: Removed column from first partition\n')
else
    %fprintf('BuildY: Removed column from second partition\n')
end

cols = m+1;
rows = m+n-t+1;

% Get the binomial coefficients corresponding to v
bi_nt = zeros(1,n-t+1);
for i = 0:1:n-t
    bi_nt(i+1) = nchoosek(n-t,i);
end

% Get the binomial coefficients corresponding to u
bi_mt = zeros(1,m-t+1);
for i=0:1:m-t
    bi_mt(i+1) = nchoosek(m-t,i);
end

% concatenate the two sets of binomial coefficients
bi = [bi_nt bi_mt];

% % Insert a zero into the x  vector at position 'q'

% get the first part of the x vector
xa = x(1:q-1) ;
% get the second part of the x vector
xb = x(q:end) ;
x = [xa; 0 ;xb] ;% Insert zero into vector

x_bi = x.*bi';


temp_matrix_1 = (zeros(rows,cols));

for j = 0:1:m
    for i = j:1:j+(n-t)
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if (i-j) == q-1
        else
            temp_matrix_1(i+1,j+1) = x_bi(i-j+1) ;
        end
    end
end

Qm = zeros(1,m+1);
for i = 0:1:m
    Qm(i+1) = nchoosek(m,i) .* (theta^i) ;
end

Qm = diag(Qm);

temp_matrix_1 = temp_matrix_1 * Qm;

temp_matrix_2 = (zeros(rows,cols));
% for each column
for j = 0:1:m-1
    % for each row
    for i = j:1:j+(m-t)
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if (i-j + (n-t+1)) == q - 1
        else
            
            temp_matrix_2(i+1,j+1) = x_bi(i-j+1 + (n-t+1));
        end
    end
end

Qn = zeros(1,n+1);
for i = 0:1:n
    Qn(i+1) = nchoosek(n,i) .* (theta^i) ;
end
Qn_a = diag([Qn 0]);
Qn_b = diag([0 Qn]);

temp_matrix_2 = temp_matrix_2 * Qn_a;



temp_matrix_3 = (zeros(rows,cols));
% for each column

for j = 1:1:m
    for i = j-1:1:j+(m-t)-1
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if i-j + (n-t+1) +1 == q - 1
            %temp_matrix_3(i+1,j+1) = 9999 ;
        else
            temp_matrix_3(i+1,j+1) =  x_bi(i-j+1 + (n-t+1) + 1) ;
        end
    end
    
end

temp_matrix_3 = temp_matrix_3 * Qn_b;

D = zeros(1,m+n-t+1);
for i = 0:1:m+n-t
    D(i+1) = 1./nchoosek(m+n-t,i);
end
D = diag(D);



Y = D*(temp_matrix_1 - (alpha.*m.*ratio.*temp_matrix_2) + (alpha.*m.*ratio.*temp_matrix_3));

end