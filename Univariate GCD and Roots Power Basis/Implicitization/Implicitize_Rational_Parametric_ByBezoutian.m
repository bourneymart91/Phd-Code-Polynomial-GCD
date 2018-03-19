function [] = Implicitize_Rational_Parametric_ByBezoutian()
% This work is based on work by Monacha - Algorithms for Intersections 1.

% Algorithm 1 for computing the entries of M (Manocha Demmel 1994)
% 1. Set index to zero
% 2. for (

x_t = [4 -3 6 6];
y_t = [1 6 0 -4];
w_t = [1 6 6 1];

m = 3;

M1 = getM(w_t,y_t,m)
M2 = getM(x_t,w_t,m)
M3 = getM(x_t,y_t,m)

X = sym('X')
Y = sym('Y')

M = X.*M1 + Y.*M2 + M3
det(M)

Implicitize_Rational_Parametric(x_t',y_t',w_t')




end

function M = getM(f,g,m)

% 1. Index = 0.
index = 0;

mypoly = []

% 2. for (0 <= j <=i) do
for i = 0:1:m
    
    % (a) for (0 <= j <= i) do
    for j = 0:1:i
        
        % i. Monomial.coeff = F[i]*G[j]
        Monomial.Coeff = f(i+1) * g(j+1);
        
        % ii. bound = i - j + 1.
        bound = i - j - 1;
        
        % iii. for ( 0 <= k <= bound) do
        for k = 0:1:bound
            
            % A. Monomials.s = i - 1 - k.
            Monomial.s = i - 1 - k;
            
            % B. Monomial.t = k + j.
            Monomial.t = k + j;
            
            % C. Poly[index ++] = Monomial.
            mypoly(index+1,1,1) = Monomial.Coeff;
            mypoly(index+1,1,2) = Monomial.s;
            mypoly(index+1,1,3) = Monomial.t;
            
            index = index + 1;
        end
    end
    
    
    
    % (b) for (i+1 <= j <= m) do
    for j = i+1:1:m
        
        %i. monomial.coeff = -F[i] * G[j]
        Monomial.Coeff = -f(i+1) * g(j+1);
        
        %ii. bound = j - i - 1
        bound = j - i - 1;
        
        %iii. for (0 <= k <= bound) do
        for k = 0:1:bound
            
            % A. Monomial.s = j - 1 - k.
            Monomial.s = j - 1 - k;
            
            % B. Monomial.t = k + i.
            Monomial.t = k+i
            
            % C. Poly[index ++] = Monomial
            mypoly(index+1,1,1) = Monomial.Coeff;
            mypoly(index+1,1,2) = Monomial.s;
            mypoly(index+1,1,3) = Monomial.t;
            
            index = index + 1;
        end
    end
    
    
end

display(mypoly)

P = zeros(m,m)

% 3. for (0 <= i <= index) do
for i = 0:1:index -1
    
    % (a) Monomial = Poly[i]
    Monomial.Coeff = mypoly(i+1,1,1)
    Monomial.s = mypoly(i+1,1,2)
    Monomial.t = mypoly(i+1,1,3)
    
    % (b) j = Monomial.s
    j = Monomial.s
    
    % (c) k = Monomial.t
    k = Monomial.t
    
    % (d) P[j][k] = P[j][k] + Monomial.coeff
    P(j+1,k+1) = P(j+1,k+1) + Monomial.Coeff
end

M = P

end