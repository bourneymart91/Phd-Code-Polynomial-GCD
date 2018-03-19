
function []=Euclid(ex_num, ec, tolerence)

% This function implements the Euclidean algorithm for Bernstein
% basis polynomials. The algorithm is described in the paper 

% Algorithm 812:BPOLY: An object-oriented library of numerical  
% algorithms for polynomials in Bernstein form

% Y. Tsai and R. Farouki, ACM Transactions on Mathematical Software, 
% volume 27, number 2, June 2001, pp. 267-296.

% n          :  An integer that defines the pair of polynomials in the 
%               database in the programme ex.m. 

% ec         :  The ratio (noise level)/(signal level) in the 
%               componentwise sense. 

% tolerance  :  The tolerance for the divisions in Euclids's algorithm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b,exactGCD,exactdeg] = ex(ex_num);
% On exit, the matrices a and b define the polynomials whose GCD is 
% to be computed. The matrix exactGCD is the GCD of the polynomials, 
% and d is the degree of the GCD of the polynomials.

% Write out the roots and their multiplicities. First write the 
% coefficients of a, and then write the coefficients of b.
disp (' ');
disp ('--------------------------------------------------');
disp (' ');
disp ('    exact root      multiplicity');
disp ('--------------------------------');
for i=1:size(a,1)
    fprintf('% 16.8e %10.0f\n',a(i,:));
end
disp (' ');

% Write the coefficients of b.
disp (' ');
disp ('    exact root      multiplicity');
disp ('--------------------------------');
for i=1:size(b,1)
    fprintf('% 16.8e %10.0f\n',b(i,:));
end
disp ('--------------------------------');
disp (' ');

disp ('The degree of the exact GCD is');
disp (exactdeg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Compute the coefficients of the polynomials from the arrays a and b.
% Note that the coefficients of the polynomils in the scaled Bernstein 
% basis are returned. 
F=B_poly(a);
G=B_poly(b);

% Calculate the coefficients of the Bernstein basis forms of
% the polynomials.
fx = sub(F);
gx = sub(G);

% Add noise to the coefficients of the Bernstein basis polynomials.
fx_noisy = noise(fx, ec, 23);
gx_noisy = noise(gx, ec, 145);

% Normalise the coefficients of f and g.
fx_noisy=fx_noisy ./ Norm(fx_noisy);
gx_noisy=gx_noisy ./ Norm(gx_noisy);

% Copy f and g into two arrays that will be required later. 
B = fx_noisy;
r = gx_noisy;

% Set the tolerence and the flag cond to detemine if the loop 
% for Euclid's algorithm should terminate.
t=tolerence;
cond=0;

% Set the counter for the number of iterations. 
Iteration_number=0;

% Implement the division sequence.
while cond==0
    
    Iteration_number=Iteration_number+1;
    
    % Perform the divisions in Euclid's algorithm.
    [B, r]=division(B,r);
    
    % The results of these two divisions determine if the 
    % algorithm is to be terminated.
    [~,r1] = division(fx_noisy, B);
    [~,r2] = division(gx_noisy, B);
    
    % Check if the algorithm is to be terminated.
    if (Norm(r1)<t) && (Norm(r2)<t)
       cond=1; 
    end 
    
end  

% Normalise the computed GCD by its 2-norm and calculate its degree.
CGCD = B / norm(B);
computeddeg = length(CGCD) - 1;

disp (' ');
disp ('The degree of the computed GCD is');
disp (computeddeg);

disp (' ')
disp ('Number of iterations:')
disp (Iteration_number);

% Print out the error in the computed GCD if the degree of the
% exact and computed GCDs are equal.
if (exactdeg == computeddeg)   

    % Calculate the coefficients of the exact GCD, and then
    % normalise the GCD by the 2-norm of their coefficients.
    C=B_poly(exactGCD);
    exGCD=sub(C);
    exGCD=exGCD/norm(exGCD);
    
    % Compute the error between the exact and computed GCDs. It is 
    % first necessary to check that the vectors of coefficients are 
    % parallel, and not anti-parallel.
    % Use the element of maximum magnitude of the exact GCD as the test.
    [~,loc]=max(abs(exGCD)); 
    if (exGCD(loc)*CGCD(loc)) < 0  % the vectors are anti-parallel
       CGCD=-CGCD;  
    end
    
    % Calculate the error.
    error = norm(CGCD-exGCD);    
    disp(' ')
    disp('The relative error in the computed GCD is');
    fprintf('% 16.8e\n ',error);
    disp(' ')
    
end

