function []  = Test

syms b_00 b_01 b_02 b_03
syms b_10 b_11 b_12 b_13
syms b_20 b_21 b_22 b_23
syms b_30 b_31 b_32 b_33

degree = 3;
switch degree
    case 1
        
        fxy = [...
            b_00 b_01
            b_10 0
            ];
    case 2
        fxy = [...
            b_00 b_01   b_02
            b_10 b_11   0
            b_20 0      0
            ];
    case 3
        fxy = [...
            b_00    b_01   b_02    b_03
            b_10    b_11   b_12    0
            b_20    b_21   0       0
            b_30    0      0       0
            ];
end

% Get the total degree
[nRows,nCols] = size(fxy);

tdeg = nRows -1;
n = tdeg;

x = sym('x');
y = sym('y');

sum = 0;

for i = 0:1:tdeg
    for j = 0:1:tdeg
        
        if i+j<=n
            
            sum = sum + ...
                (...
                fxy(i+1,j+1) * Trinomial(n,i,j) *...
                x^(i) * y^(j) * (1-x-y)^(n-i-j)...
                );
        end
        
    end
end

my_exp = expand(sum)
collect(my_exp,[x,y])

matrix_2 = [...
     1  0  0 0 0 0 
    -2  1  0 0 0 0
    -2  0  1 0 0 0
     1 -1  0 1 0 0
     2 -1 -1 0 1 0
     1 -1  0 0 0 1
    ];

trinoms_2 = BuildQ(2);

matrix_2_bi = matrix_2 * trinoms_2;

pinv(matrix_2_bi)


matrix_3 = [...
     1  0  0  0  0  0 0 0 0 0
    -3  3  0  0  0  0 0 0 0 0
    -3  0  3  0  0  0 0 0 0 0
     3 -6  0  3  0  0 0 0 0 0
     6 -6 -6  0  6  0 0 0 0 0
     3  0 -6  0  0  3 0 0 0 0
    -1  3  0 -3  0  0 1 0 0 0
    -3  6  3 -3 -6  0 0 3 0 0
    -3  3  6  0 -6 -3 0 0 3 0
    -1  0  3  0  0 -3 0 0 0 1
 ];





% 
% % Get set of basis elements
% nBasisElements = nchoosek(n+2,2);
% 
% %mat_basis = zeros(n+1,n+1);
% 
% for i = 0:1:tdeg
%     for j = 0:1:tdeg
%         
%         if i+j <=n
%             xx = (factorial(n)./ (factorial(i)*factorial(j)*factorial(n-i-j))) ...
%                 * x^(i) * y^(j) * (1-x-y)^(n-i-j);
%             
%             display(xx);
%         end
%         
%     end
% end


end

function Q = BuildQ(m)

Tri_vec = GetTrinomials(m);

Q = diag(Tri_vec);


end










