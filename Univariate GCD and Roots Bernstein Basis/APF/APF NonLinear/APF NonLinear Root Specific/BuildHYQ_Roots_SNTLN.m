function H1Y1Q1 = BuildHYQ_Roots(dx,m,theta)
% Used in APF_Roots
% See Report (APF - Roots - Rearrangement of the Matrix Vector Product)
%
% Build HYQ
% where H is a diagonal matrix of binomials corresponding to f and g
% where Y is a Toeplitz matrix in terms of d
% where Q is a diagonal matrix of binomials corresponding to u and v
%
% % Inputs
%
% dx : Coefficients of polynomial d(x)
%
% m : Total degree of the polynomial f(x)
%
% theta : Optimal value of \theta
%
% % Outputs
%
% H1Y1Q1 :


H1Y1Q1 = BuildHYQ_Roots2(dx,m,theta);


end


function H1Y1Q1 = BuildHYQ_Roots1(dx,m,theta)

% Used in APF_Roots
% See Report (APF - Roots - Rearrangement of the Matrix Vector Product)

% Build HYQ
% where H is a diagonal matrix of binomials corresponding to f and g
% where Y is a Toeplitz matrix in terms of d
% where Q is a diagonal matrix of binomials corresponding to u and v

% HYQ*[u;v] is equivalent to HCG*d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                         Inputs


%   dx - coefficients of the GCD

%   m - degree of polynomial m

%   theta - optimal value of \theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the degree of polynomial dx
t = length(dx)-1;

% Initialise an empty matrix.
H1Y1Q1_Partition = zeros(m+1,m-t+1);

% for each column 0:1:m-t
for j = 0:1:m-t
    for i = j:1:j+t
        %   HY_Partition1(i+1,j+1) = dx(i-j+1) .* theta^(i-j) .* nchoosek(t,i-j) ./ nchoosek(m,i);
        H1Y1Q1_Partition(i+1,j+1) = dx(i-j+1) .* theta^(i) .* nchoosek(t,i-j) .*nchoosek(m-t,j) ./ nchoosek(m,i);
    end
end

H1Y1Q1 = H1Y1Q1_Partition;

end

function H1Y1Q1 = BuildHYQ_Roots2(dx,m,theta)
% Used in APF_Roots
% See Report (APF - Roots - Rearrangement of the Matrix Vector Product)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                         Inputs


%   dx

%   m

%   theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Get the degree of polynomial dx
t = length(dx)-1;

% Build Matrix H
H1 = zeros(m+1,1);

for i=0:1:m
    H1(i+1) = nchoosek(m,i);
end

H1 = diag(1./H1);

% Build Matrix Y1
Y1 = zeros(m+1,m-t+1);
for j = 0:1:m-t
    for i = j:1:j+t
        %   HY_Partition1(i+1,j+1) = dx(i-j+1) .* theta^(i-j) .* nchoosek(t,i-j) ./ nchoosek(m,i);
        Y1(i+1,j+1) = dx(i-j+1) .* nchoosek(t,i-j) ;
    end
end


% Build thetas matrix
theta_vector = zeros(m+1,1);
for i = 0:1:m
    theta_vector(i+1) = theta^(i);
end
theta_matrix = diag(theta_vector);

% build Q consisting of binomial coefficients corresponding to u
Q1 = zeros(m-t+1,1);

for i = 0:1:m-t
    Q1(i+1) = nchoosek(m-t,i);
end

Q1 = diag(Q1);

H1Y1Q1 = H1*theta_matrix*Y1*Q1;


end