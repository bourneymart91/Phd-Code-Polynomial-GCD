% function y = LSE(E,f,C,g)
% % This function uses the QR decomposition to solve the LSE problem
% % minimise ||Ey-f||  subject to  Cy=g
% % The output is the vector y.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     [m1,~] = size(C);
%     [Q,R] = qr(C'); 
%     R1 = R(1:m1,:);
%     w1 = R1'\g;
%     %Q1 = Q(:,1:m1);
%     Q2 = Q(:,m1+1:end);
% 
%     w2 = Q2'*f;
% 
%     y = Q*[w1;w2];
% 
% 
% 
% end


function x = LSE(A,b,B,d)

% This function solves the least square equality problem.

%   minmize ||Ax-b||_{2}, subject to the constraint Bx=d

% using Algorithm 12.1.2, page 586, in Matrix Computations 
% by G. Golub and C. Van Loan, 1996.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p, n]=size(B);

% Step 1:   The QR factorization of B' (the transpose of B)
[Q,R]=qr(B');
R1=R(1:p,1:p);

% Step 2:    Solve R1'y=d for y
y=R1'\d;

% Step 3:    
% Let A=AQ, then find z so ||A(:,p+1:n)z-(b-A(:,1:p)y)|| is minimised
A = A*Q;
A1 = A(:,1:p);
A2 = A(:,p+1:n);
z = pinv(A2)*(b-A1*y);

% Step 4: x=Q(:,1:p)y+Q(:,p+1:n)z
x = Q*[y;z];

end


