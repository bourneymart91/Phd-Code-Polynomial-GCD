function [] = test_APF()


m = 5;
n = 4;
t = 2;
m_t = m-t;
n_t = n-t;

% create symbolic vector of coefficients of u
u = sym('U',[m_t+1,1]);
d = sym('D',[t+1,1]);

d_binoms = sym(zeros(t+1,1))

for i = 0:1:t
   d_binoms(i+1) =nchoosek(t,i); 
end

% C*d should be equal to Yu
d_binoms = zeros(t+1,1)
for i = 0:1:t
   d_binoms(i+1) = nchoosek(t,i) ;
end
d_binoms = diag(d_binoms)

C = BuildC(u,m_t,t,m,n);
Cd = expand(C*d_binoms*d)

Y = BuildY(d_binoms*d,m,t)

u_binoms = zeros(m-t+1,1);

for i = 0:1:m_t
   u_binoms(i+1) = nchoosek(m_t,i); 
end

u_binoms = diag(u_binoms);


YQU = expand(Y*(u_binoms*u))
Cd
YQU./Cd
end


function [C] =  BuildC(u,m_t,t,m,n)

C1 = sym(zeros(m+1,t+1));


% % Build the matrix C1 without binomial coefficients
% for each column 1,...,t+1
for j = 0:1:t
    % for each coefficient of u u_{0},...,u_{m-t}
    for i = j:1:j+m_t
        C1(i+1,j+1) = u(i-j+1) .* nchoosek(m-t,i-j) ./ nchoosek(m,i);
    end
end

[rows,cols] = size(C1)

% C2 Part 1 is given by C1 with the first row removed.
A1 = [zeros(rows-1,1) eye(rows-1)];
C2_Part1 = A1 * C1;
A2 = [eye(rows-1) zeros(rows-1,1)];
C2_Part2 = A2 * C1;

C2 = m.* (C2_Part1 - C2_Part2)

% % Concatenate C1 and C2

C = [C1 ; C2]

end



function [Y] = BuildY(d_bi,m,t)

n = m-1;

% % Build matrix Y1 Excluding the binomials of H^-1 and Q

% Initialise matrix Y1
Y1 = sym(zeros(m+1,m-t+1));

% for each column
for j = 0:1:m-t
    % for each coefficient
    for i = j:1:j+t
        Y1(i+1,j+1) = d_bi(i-j+1) ./ nchoosek(m,i) ;
    end
end

% % Build matrix Y2 Excluding the binomials of H^-1 and Q
Y1
% get number of rows in Y1
[rows,cols] = size(Y1);

% Let A1 be the identity matrix with a zero column such that A*Y gives the
% matrix Y with the first row removed.


A1 = [zeros(rows-1,1) eye(rows-1)];
Y2_Part1 = m.*A1 * Y1

A2 = [eye(rows-1) zeros(rows-1,1)];
Y2_Part2 = m.*A2 * Y1

Y2 =  (Y2_Part1 - Y2_Part2)

% % Concatenate Y1 and Y2

Y = [Y1; Y2];

end