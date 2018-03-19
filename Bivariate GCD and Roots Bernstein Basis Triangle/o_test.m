% Test QR decompositon
% Test computing sk from sk-1

% Add subfolders
restoredefaultpath

addpath (...
        'Examples',...
        'Formatting',...
        'GetCofactors',...
        'GetGCDCoefficients',...
        'GetGCDDegree',...
        'Low Rank Approximation',...
        'Preprocessing'...
        );
    
    
% Get Random polynomials f(x,y) and g(x,y)
ex_num = '1';
[fxy,gxy,d,u,v] = Examples_GCD(ex_num);

m = GetDegree_Bivariate(fxy);
n = GetDegree_Bivariate(gxy);

for k = 1 : 1 : min(m,n)

    S{k} = BuildDTQ_2Polys(fxy,gxy,k);
    
    if k >= 2
        S_calc{k} = BuildDTQ_from_prev(S{k-1},m,n,k);
        

    else
        S_calc{1} = S{1};
        
    end
    
    vSingularValues(k) = min(svd(S{k}));
    vSingularValues_calc(k) = min(svd(S_calc{k}));
    
    [Q{k},R{k}] = qr(S{k});
    
    [nRows_Rk,nCols_Rk] = size(R{k});
    R1{k} = R{k}(1:nCols_Rk,1:nCols_Rk);
    diag_R1{k} = diag(R1{k});
    display(diag_R1{k})
    display('end\n');
    
end

% Plot the unordered diagonals of R
figure()
hold on 
for k = 1:1:min(m,n)
    nEntries = length(diag_R1{k});
    plot((1:1:nEntries),log10(diag_R1{k})','-s') 
end
hold off

figure()
hold on
for k = 1:1:min(m,n)
    nEntries = length(diag_R1{k});
    
    plot((1:1:nEntries),sort(log10(diag_R1{k}))','-s') 
end
hold off


figure()
plot(log10(vSingularValues),'-s')
hold on
plot(log10(vSingularValues_calc),'-s')
hold off

% Get singular values of

