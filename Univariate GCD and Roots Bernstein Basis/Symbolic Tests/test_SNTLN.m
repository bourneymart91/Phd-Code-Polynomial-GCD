function [] = test_SNTLN()



% Set degree of input polynomial f
m = 4;

% set degree of derivative polynomial g
n = m-1;

% set degree of the subresultant S_{t}
t = 1;

% Set the index of the column to be removed starting at 0
q = 1;

% create symbolic vector of coefficients
f = sym('A',[m+1,1])

% create symbolic theta and alpha
th = sym('th');
al = sym('al');

% build the symbolic sylvester subresultant by old method
S_Brnstn_old = buildSymbolicSylvesterMatrix_Bernstein_old(f,m,n,t,th,al)
% Build the symbolic sylvester subresultant by new method
S_Brnstn_new = buildSymblcSylvstrMatrix_Brnstn_new(f,m,n,t,th,al)

% get the column to be removed
if q >= 0
    S_remvd_col = S_Brnstn_old(:,q+1);
end
% Build the symbolic sylvester matrix without binomial coefficients
% S_pwr = buildSymblcSylvstrMtrx_Pwr(f,m,n,t)

% build the symbolic vector x such that S_{k} x_{k} = 0
x = sym('X',[m+n-2*t+2,1])

% Remove the column from the Sylvester matrix
% Ak_pwr = S_pwr;
% if q >=0
%     Ak_pwr(:,q+1) = [];
% end

Ak_Brnstn = S_Brnstn_new
if q>=0
    ck_Brnstn = Ak_Brnstn(:,q+1);
    Ak_Brnstn(:,q+1) = []
end

% remove corresponding value from the x vector
fprintf('x vector of Ak \n')
x_removed_col = x;
if q>=0
    x_removed_col(q+1) = [];
end


% fprintf('Build Matrix Y without binomials\n')
% % get result S
% % Build matrix Y such that Y*f = S*x
% Y = buildY(f,m,n,t,x,q)
% 
% % Build ck by method of sylvester matrix times x
% ck = Ak_pwr * x_removed_col
% % Build ck by rearrangment method
% ck2 = Y * f
% % Get difference between the two vectors - Should be zero vector
% ck - ck2
% % j ranges from 0,...,m+n-2t+1
% % numer of columns = m+n-2t+2

fprintf('Build Matrix Y with binomials\n')
% get result S
% Build matrix Y such that Y*f = S*x
Y = buildY_binoms(m,n,t,q,x,al,th)

% Build ck by method of sylvester matrix times x
Ax_Brnstn = Ak_Brnstn * x_removed_col ;
Ax_Brnstn = expand(Ax_Brnstn)

% Build ck by rearrangment method
Ax2_Brnstn = Y * f ;
Ax2_Brnstn = expand(Ax2_Brnstn)


Ax_Brnstn - Ax2_Brnstn

% Get difference between the two vectors - Should be zero vector

% j ranges from 0,...,m+n-2t+1
% numer of columns = m+n-2t+2


if q < 0
else
    fprintf('Build P using old method\n')
    fprintf('Pf = ck \n')
    P = buildP_old(f,m,n,t,q,al,th);
    
    ans1 = P*f
    
    
    fprintf('Build P using new method\n')
    fprintf('Pf = ck \n')
    P = buildP_new(f,m,n,t,q,al,th);
    ans2 = P*f
    
    fprintf('ck is given by')
    S_remvd_col = expand(S_Brnstn_old(:,q+1))
    expand(ck_Brnstn)
end
end


function [Y] = buildY(f,m,n,t,x,q)

cols = m+1;
rows = m+n-t+1;

%temp_matrix_1 = zeros(m+n-t+1, m+1)
for j = 0:1:m
    for i = j:1:j+(n-t)
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if (i-j) == q
        else
            temp_matrix_1(i+1,j+1) = x(i-j+1);
        end
    end
end
temp_matrix_1

temp_matrix_2 = sym(zeros(rows,cols));
% for each column
for j = 0:1:m-1
    % for each row
    for i = j:1:j+(m-t)
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if (i-j + (n-t+1)) == q
        else
            temp_matrix_2(i+1,j+1) = m*x(i-j+1 + (n-t+1));
            
        end
    end
end

temp_matrix_2 = temp_matrix_2

temp_matrix_3 = sym(zeros(rows,cols));
% for each column
for j = 1:1:m
    for i = j-1:1:j+(m-t)-1
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if i-j + (n-t+1) == q
        else
            temp_matrix_3(i+1,j+1) = m*x(i-j+1 + (n-t+1) + 1);
        end
    end
    
end
temp_matrix_3

Y = temp_matrix_1 - temp_matrix_2 + temp_matrix_3;

end

function [Y] = buildY_binoms(m,n,t,q,x,al,th)

cols = m+1;
rows = m+n-t+1;

bi_nt = zeros(1,n-t+1);
for i = 0:1:n-t
    bi_nt(i+1) = nchoosek(n-t,i);
end
bi_mt = zeros(1,m-t+1);
for i=0:1:m-t
    bi_mt(i+1) = nchoosek(m-t,i);
end

bi = [bi_nt bi_mt];
x_bi = x.*bi';


%temp_matrix_1 = zeros(m+n-t+1, m+1)
for j = 0:1:m
    for i = j:1:j+(n-t)
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if (i-j) == q
        else
            %temp_matrix_1(i+1,j+1) = x(i-j+1) .*nchoosek(n-t,i-j);
            temp_matrix_1(i+1,j+1) = x_bi(i-j+1) ;
        end
    end
end

for i = 0:1:m
    Qm(i+1) = nchoosek(m,i) .* th^i ;
end

Qm = diag(Qm);

temp_matrix_1 = temp_matrix_1 * Qm

temp_matrix_2 = sym(zeros(rows,cols));
% for each column
for j = 0:1:m-1
    % for each row
    for i = j:1:j+(m-t)
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if (i-j + (n-t+1)) == q
        else
            %temp_matrix_2(i+1,j+1) = m*x(i-j+1 + (n-t+1)) .* nchoosek(m-t,i-j);
            temp_matrix_2(i+1,j+1) = al.* m .* x_bi(i-j+1 + (n-t+1));
        end
    end
end

for i = 0:1:n
    Qn(i+1) = nchoosek(n,i) .* th^i ;
end
Qn_a = diag([Qn 0]);
Qn_b = diag([0 Qn]);

temp_matrix_2 = temp_matrix_2 * Qn_a



temp_matrix_3 = sym(zeros(rows,cols));
% for each column

for j = 1:1:m
    for i = j-1:1:j+(m-t)-1
        % if i-j gives the index of the removed column then leave matrix
        % entry as zero, otherwise...
        if i-j + (n-t+1) == q
        else
            %temp_matrix_3(i+1,j+1) = m*x(i-j+1 + (n-t+1) + 1) .* nchoosek(m-t,i+1-j);
            temp_matrix_3(i+1,j+1) = al .* m .* x_bi(i-j+1 + (n-t+1) + 1) ;
        end
    end
    
end

temp_matrix_3 = temp_matrix_3 * Qn_b


for i = 0:1:m+n-t
    D(i+1) = 1./nchoosek(m+n-t,i);
end
D = diag(D);

Y = D* [temp_matrix_1 - temp_matrix_2 + temp_matrix_3];

end

function [S] = buildSymblcSylvstrMtrx_Pwr(f,m,n,t)

C1 = buildSymblcSylvstrMtrxPartition_Pwr(f,m,n,t)
C2 = buildSymblcSylvstrMtrxPartition2_Pwr(f,m,n,t)
S = [C1 C2];
end

function [C1] = buildSymblcSylvstrMtrxPartition_Pwr(f,m,n,t)
for j = 0:1:n-t
    % for each coefficient a in f
    for i = j:1:j+m
        C1(i+1,j+1) = f(i-j+1).* (th^(i-j));
    end
end
end
function [C2] = buildSymblcSylvstrMtrxPartition2_Pwr(f,m,n,t,th,al)
for j = 0:1:m-t
    % for each coefficient a in f
    for i = j:1:j+n
        C2(i+1,j+1) = al.* m*(f(i-j+2) - f(i-j+1)) .* (th^(i-j));
    end
end
end

function [S] = buildSymblcSylvstrMatrix_Brnstn_new(f,m,n,t,th,al)
% uses the sylvester rearrangement
C1 = buildSymbolicSylvesterPartition_Bernstein_new(f,m,n,t,th,al);
C2 = buildSymbolicSylvesterPartition2_Bernstein_new(f,m,n,t,th,al);

S = [C1 C2];
end

function [S] = buildSymbolicSylvesterMatrix_Bernstein_old(f,m,n,t,th,al)
% uses the old format D^{-1}TQ

C1 = buildSymbolicSylvesterPartition_Bernstein_old(f,m,n,t,th,al);
C2 = buildSymbolicSylvesterPartition2_Bernstein_old(f,m,n,t,th,al);

S = [C1 C2];
end

function [C] = buildSymbolicSylvesterPartition_Bernstein_old(f,m,n,t,th,al)

%C = zeros(m+n-t+1,n-t+1);

% for each column
for j = 0:1:n-t
    % for each coefficient a in f
    for i = j:1:j+m
        C(i+1,j+1) = f(i-j+1) .*(th^(i-j)) .* nchoosek(m,i-j) .* nchoosek(n-t,j) ./ nchoosek(m+n-t,i);
    end
end
end

function [C] = buildSymbolicSylvesterPartition2_Bernstein_old(f,m,n,t,th,al)

for j = 0:1:m-t
    % for each coefficient a in f
    for i = j:1:j+n
        C(i+1,j+1) = al.* m*(f(i-j+2) - f(i-j+1)) .* (th^(i-j)) .* nchoosek(n,i-j) .* nchoosek(m-t,j) ./ nchoosek(m+n-t,i);
    end
end


end

function [C] = buildSymbolicSylvesterPartition_Bernstein_new(f,m,n,t,th,al)
% f - coefficients of polynomial f

% m - degree of polynomial f

% n - degree of polynomial g

% th - optimal value of theta

% al - optimal value of alpha

% for each column
for j = 0:1:n-t
    % for each coefficient a in f
    for i = j:1:j+m
        C(i+1,j+1) = f(i-j+1).*(th^(i-j)) .* nchoosek(m+n-t-i,n-t-j) .* nchoosek(i,j) ./ nchoosek(m+n-t,n-t);
    end
end


end

function [C] = buildSymbolicSylvesterPartition2_Bernstein_new(f,m,n,t,th,al)
% for each column
for j = 0:1:m-t
    % for each coefficient a in f
    for i = j:1:j+n
        C(i+1,j+1) = al.* m*(f(i-j+2) - f(i-j+1)) .* (th^(i-j)) .* nchoosek(m+n-t-i,m-t-j) .* nchoosek(i,j) ./ nchoosek(m+n-t,m-t);
    end
end

end

function [P] =  buildP_old(f,m,n,t,q,alpha,theta)

if q <= n-t
    P = buildP_LHS_old(f,m,n,t,q,alpha,theta);
elseif q<m+n-(2*t)+2
    P = buildP_RHS_old(f,m,n,t,q,alpha,theta);
end

end

function [P] = buildP_new(f,m,n,t,j,alpha,theta)


if j <= n-t
    P = buildP_LHS_new(f,m,n,t,j,alpha,theta);
elseif j<m+n-2*t+2
    P = buildP_RHS_new(f,m,n,t,j,alpha,theta);
end




end


function [P] = buildP_LHS_new(f,m,n,t,q,alpha,th)


Z1 = zeros(q,m+1);

x = (m+n-t+1) - q - (m+1);
Z2 = zeros(x,m+1);

G = sym(zeros(m+1,m+1));
for i = 0:1:m
    G(i+1,i+1) = 1 .* th^i .* nchoosek(i+q,q) .* nchoosek(m+n-t-(i+q),m-i) ./ nchoosek(m+n-t,n-t);
end
size(G)
P = [Z1 ; G; Z2];



end

function [P] = buildP_LHS_old(f,m,n,t,q,alpha,th)

Z1 = zeros(q,m+1);

x = (m+n-t+1) - q - (m+1);
Z2 = sym(zeros(x,m+1));

G = sym(zeros(m+1,m+1));
for i = 0:1:m
    G(i+1,i+1) = nchoosek(m,i) .*(th^i) .* nchoosek(n-t,q) ./nchoosek(m+n-t,q+i);
end
size(G)
P = [Z1 ; G; Z2];



end

function [P] = buildP_RHS_new(f,m,n,t,q,al,th)

j = q - (n-t+1);

Z1 = zeros(j,m+1);

y = (m+n-t+1)- j - m;

Z2 = zeros(y,m+1);

G = sym(zeros(m,m+1));

% for each of the coefficients of c_{q}
for i = 0:1:m-1
    G(i+1,i+1) = - al.* m .* th^i .*  nchoosek(i+j,j) .* nchoosek(m+n-t-(i+j),n-i) ./ nchoosek(m+n-t,m-t);
    G(i+1,i+2) =   al.* m .* th^i .*  nchoosek(i+j,j) .* nchoosek(m+n-t-(i+j),n-i) ./ nchoosek(m+n-t,m-t);
    
    
end
size(G)
P = [Z1; G; Z2];
end


function [P] = buildP_RHS_old(f,m,n,t,q,al,th)

j = q-(n-t+1);

Z1 = zeros(q-(n-t+1),m+1);

y = (m+n-t+1)-(q-(n-t+1)) - m;

Z2 = zeros(y,m+1);

G = sym(zeros(m,m+1));

for i = 0:1:m-1
    G(i+1,i+1) = - al.* m .*(th^i) .* nchoosek(n,i) .*nchoosek(m-t,j) ./ nchoosek(m+n-t,j+i);
    G(i+1,i+2) =   al.* m .*(th^i) .* nchoosek(n,i) .*nchoosek(m-t,j) ./ nchoosek(m+n-t,j+i);
end

P = [Z1; G; Z2];
end