r
function []=o(fx_noisy, gx_noisy)

% Normalise the coefficients of f and g by the geometric means of
% their coefficients. Repeat this procedure for noisef and noiseg.
GM_fx = GM(fx_noisy);
GM_gx = GM(gx_noisy);

fx_noisy = fx_noisy/GM_fx;
gx_noisy = gx_noisy/GM_gx;



% Degree elevate the polynomial of lower degree in preparation for the
% construction of the Bezout matrix.
m = length(fx_noisy) - 1;
n = length(gx_noisy) - 1;

if m>n
    gx_noisy = deg_elevation(gx_noisy,n,m);
    noise_gx = deg_elevation(noise_gx,n,m);
elseif m<n
    fx_noisy = deg_elevation(fx_noisy,m,n);
    noise_fx = deg_elevation(noise_fx,m,n);
end

% Construct the Bezout matrix B for theta=1. Do this for (f,g)
% and (noisef,noiseg).
B_fg = bezout(fx_noisy, gx_noisy);
noise_Bfg = bezout(noise_fx, noise_gx);

% Calculate the 2-norms of B and noiseB, and write out their values.
normB = norm(B_fg);
normnoiseB = norm(noise_Bfg);

disp('----------------------------------------------------');
disp(' ');
disp('The 2-norm of the Bezout matrix of the inexact polynomials');
disp('for theta=1 is:');
fprintf('% 16.8e \n',normB);
disp(' ');
disp('The 2-norm of the Bezout matrix of the perturbation polynomials');
disp('for theta=1 is:');
fprintf('% 16.8e \n',normnoiseB);
disp(' ');
disp('----------------------------------------------------');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now consider the case theta is not equal to one.

% Calculate the optimal value of theta for transformation from the
% Bernstein basis to the modified Bernstein basis.
theta = optimal(B_fg);
disp(' ');
disp('The optimal value of theta is');
fprintf('%8.5e \n ',theta);
disp(' ');

% Use the optimal value of theta to calculate the Bezout matrix B2
% with theta~=1. Repeat this for the Bezout matrix B2 of the
% perturbation polynomials noisef and noiseg.
B2 = bezout_theta(B_fg,theta);
noiseB2=bezout_theta(noise_Bfg,theta);

% Calculate the 2-norms of B2 and noiseB2, and write out their values.
normB2=norm(B2);
normnoiseB2=norm(noiseB2);

disp('----------------------------------------------------');
disp(' ');
disp('The 2-norm of the Bezout matrix of the inexact polynomials');
disp('for the optimal value of theta is:');
fprintf('% 16.8e \n',normB2);
disp(' ');
disp('The 2-norm of the Bezout matrix of the perturbation polynomials');
disp('for the optimal value of theta is:');
fprintf('% 16.8e \n',normnoiseB2);
disp(' ');
disp('----------------------------------------------------');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the sum of the elements in each column of B and B2.
sumcolumnB(B_fg,1);
sumcolumnB(B2,2);

% Calculate the singular values of B, and their normalised values.
vSingularValues=svd(B_fg);
vSingularValues=log10(vSingularValues/vSingularValues(1));

% Plot the singular values. First, set the x and y axes.
n = length(vSingularValues);
x = 1 : 1 : n;
y = vSingularValues(x);

% Calculate the last non-zero singular value of B.
exact_rank=n - t_exact;
rr=y(exact_rank);   % the y-value of r for the graph.

figure(3)
plot(x,y,'b-o','MarkerSize',6,'MarkerFaceColor','b','LineWidth',1)

hold on
plot(exact_rank,rr,'rs','MarkerSize',6,'MarkerFaceColor','r')
xlabel('i')
ylabel('log_{10} \sigma_{i} / \sigma_{1}')
ts ='Bezout matrix: ';
title([ts,'\theta=1'],'FontSize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now repeat figure 1 for the singular values of B2, the Bezout matrix
% for the modified Bernstein basis.

% Prepare the data to be plotted.
vSingularValues=svd(B2);
vSingularValues=log10(vSingularValues/vSingularValues(1));
n=length(vSingularValues);

x=1:1:n;
y=vSingularValues(x);

% Calculate the last non-zero singular value of B2.
exact_rank=n-t_exact;
rr=y(exact_rank);   

figure (4)
plot(x,y,'b-o','MarkerSize',6,'MarkerFaceColor','b','LineWidth',1)

hold on
plot(exact_rank,rr,'rs','MarkerSize',6,'MarkerFaceColor','r')
xlabel('i')
ylabel('log_{10} \sigma_{i} / \sigma_{1}')
title(['Bezout matrix with \theta=',num2str(theta)],'FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the coefficients of the computed GCD from B. Reduce B to
% upper triangular form using the QR and LU decompositions. The 
% coefficients of the GCD are stored in the vectors CGCD_qr and
% CGCD_lu respectively.
[computedGCD_qr,computedGCD_lu]=gcd_B(B_fg,t_exact,1);

% Compare the theoretically exact GCD with the computed GCD.
% First, construct the vector of coefficients of the exact GCD, and
% normalise the coefficients to have unit 2-norm.
eGCD=B_poly(dx_exact);
eGCD=sub(eGCD);
eGCD=eGCD/norm(eGCD); 

% Compute the error between the exact GCD, and the GCD computed from 
% the QR and LU decompositions. It is first necessary to check that the 
% vectors of coefficients are parallel, and not anti-parallel.
% Use the element of maximum magnitude of the exact GCD as the test.
[~,loc]=max(abs(eGCD)); 
if (eGCD(loc)*computedGCD_qr(loc)) < 0  % the vectors are anti-parallel
    computedGCD_qr=-computedGCD_qr;  
end;

% Repeat for the GCD computed from the LU decomposition.
if (eGCD(loc)*computedGCD_lu(loc)) < 0  
    computedGCD_lu=-computedGCD_lu;
end;

% Compute the error between the GCD computed from the QR and LU
% decompositions.
error_qr=norm(computedGCD_qr-eGCD);
error_lu=norm(computedGCD_lu-eGCD);

disp(' ')
disp('Relative error in GCD computed from the LU decomposition');
disp('for the Bernstein basis');
fprintf('% 16.8e \n',error_lu);

disp(' ');
disp('Relative error in GCD computed from the QR decomposition');
disp('for the Bernstein basis');
fprintf('% 16.8e \n',error_qr);
disp(' ')
disp('-----------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Repeat this procedure for the modified Bernstain basis.
[computedGCD_qr, computedGCD_lu]=gcd_B(B2,t_exact,theta);

[~,loc]=max(abs(eGCD)); 
if (eGCD(loc)*computedGCD_qr(loc)) < 0  % the vectors are anti-parallel
    computedGCD_qr=-computedGCD_qr;  
end;

if (eGCD(loc)*computedGCD_lu(loc)) < 0  
    computedGCD_lu=-computedGCD_lu;
end;

error_qr=norm(computedGCD_qr-eGCD);
error_lu=norm(computedGCD_lu-eGCD);

disp(' ')
disp('Relative error in GCD computed from the LU decomposition');
disp('for the modified Bernstein basis');
fprintf('% 16.8e \n',error_lu);

disp(' ');
disp('Relative error in GCD computed from the QR decomposition');
disp('for the modified Bernstein basis');
fprintf('% 16.8e \n',error_qr);
disp(' ')

