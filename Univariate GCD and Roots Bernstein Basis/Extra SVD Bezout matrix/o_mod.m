
function [dx_qr, dx_lu] = o_mod(fx_noisy, gx_noisy, bool_preproc)

% Normalise the coefficients of f and g by the geometric means of
% their coefficients. Repeat this procedure for noisef and noiseg.
GM_fx = GM(fx_noisy);
GM_gx = GM(gx_noisy);

fx_noisy = fx_noisy/GM_fx;
gx_noisy = gx_noisy/GM_gx;


% Degree elevate the polynomial of lower degree in preparation for the
% construction of the Bezout matrix.
m = GetDegree(fx_noisy);
nSingularValues = GetDegree(gx_noisy);

m_star = max(m,nSingularValues);

if m > nSingularValues
    
    gx_noisy = deg_elevation(gx_noisy, nSingularValues, m);
    
elseif m < nSingularValues
    
    fx_noisy = deg_elevation(fx_noisy, m, nSingularValues);
    
end

% Construct the Bezout matrix B for theta=1. Do this for (f,g)
% and (noisef,noiseg).
Bez_fg = bezout(fx_noisy, gx_noisy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now consider the case theta is not equal to one.

% Calculate the optimal value of theta for transformation from the
% Bernstein basis to the modified Bernstein basis.


switch bool_preproc
    case true
        theta = optimal(Bez_fg);
    case false 
        theta = 1;
end

disp(' ');
disp('The optimal value of theta is');
fprintf('%8.5e \n ',theta);
disp(' ');

% Use the optimal value of theta to calculate the Bezout matrix B2
% with theta~=1. Repeat this for the Bezout matrix B2 of the
% perturbation polynomials noisef and noiseg.
B_preproc = bezout_theta(Bez_fg, theta);


% Now repeat figure 1 for the singular values of B2, the Bezout matrix
% for the modified Bernstein basis.

% Prepare the data to be plotted.
vSingularValues = svd(B_preproc);
vNormalisedSingularValues =  vSingularValues/vSingularValues(1);
nSingularValues = length(vNormalisedSingularValues);

xVec = 1 : 1 : nSingularValues;
y = vSingularValues(xVec);

% Calculate the last non-zero singular value of B2.
%exact_rank=n-t_exact;
%rr=y(exact_rank);   

figure (4)
plot(xVec, log10(vNormalisedSingularValues), 'b-o', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'LineWidth', 1)

hold on
%plot(exact_rank,rr,'rs','MarkerSize',6,'MarkerFaceColor','r')
xlabel('i')
ylabel('log_{10} \left( \sigma_{i} / \sigma_{1} \right)','Interpreter', 'latex')
title(['Bezout matrix with \theta=',num2str(theta)],'FontSize',12)

[~, location] = max(abs(diff(log10(vNormalisedSingularValues))));
t = m_star - location;

display(theta)

[dx_qr, dx_lu] = gcd_B(Bez_fg, t, theta);

dx_qr = dx_qr';
dx_lu = dx_lu';

dx_qr = GetWithThetas(dx_qr, theta);
dx_lu = GetWithThetas(dx_lu, theta);

end
