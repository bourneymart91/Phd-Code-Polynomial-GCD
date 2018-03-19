function B = Bezout(p,q)

% This function builds a bezoutian matrix, by the three stages defined in
% Bini's paper.

% #########################################################################

%                       Inputs.

% p - coefficients of polynomial p

% q - coefficients of polynomial q

% #########################################################################

% Get degree of polynomial p
n = length(p)-1;

% Build an empty Bezoutian matrix
B = zeros(n,n);

% Rule One.
for i = 1:1:n
    B(i,1) = (n./i) * (p(i+1)*q(1) - p(1)*q(i+1));
end

% Rule Three.
for j = 1:1:n-1
   B(n,j+1) = n./(n-j) ...
       * (p(n+1)*q(j+1) - p(j+1)*q(n+1));
end

% Rule Two.
for j=1:1:n-1
    for i=1:1:n-1
        B(i,j+1) = (n^2 ./(i*(n-j))) * ...
        (...
            p(i+1)*q(j+1) - p(j+1)*q(i+1)...
        )...
        + (j*(n-i)) ./ (i*(n-j)) * B(i+1,j);
        
    end
end

end