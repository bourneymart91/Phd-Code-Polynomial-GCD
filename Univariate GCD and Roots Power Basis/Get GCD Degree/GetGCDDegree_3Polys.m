function t = GetGCDDegree_3Polys(fx, gx, hx, limits_t, rank_range)
% GetGCDDegree(fx,gx)
%
% Get the degree of the GCD d(x) of f(x) and g(x), by Sylvester matrix method.
%
% Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficietns of polynomial g(x)
%
% hx : (Vector) Coefficients of polynomial h(x)
%
% degree_limits : [Int Int]
%
% Outputs.
%
% t : (Int) Degree of GCD of f(x) and g(x)


% Get Degree of polynomial f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% If the number of distinct roots in f(x) is one, then the degree of the
% GCD of f(x) and f'(x) = m-1 = n.

% set upper and lower bound of subresultants to be computed
lowerLimit_k = 1;
upperLimit_k = min([m,n,o]);
limits_k = [lowerLimit_k upperLimit_k];

% Get the number of subresultants which must be constructed.
nSubresultants = upperLimit_k - lowerLimit_k +1 ;

% Initialise a vector to store the minimum singular values for each
% S_{k}.
vMinimumSingularValues = zeros(1,nSubresultants);


% Initialise a vector to store the minimum distnaces for each S_{k}.
vMinimumResidual = zeros(1, nSubresultants);
vMaxDiagR1      = zeros(1, nSubresultants);
vMinDiagR1      = zeros(1, nSubresultants);
vMaxRowNormR1   = zeros(1, nSubresultants);
vMinRowNormR1   = zeros(1, nSubresultants);


% Initialise some arrays
arr_Sk = cell(1, nSubresultants);
arr_R = cell(1, nSubresultants);
arr_R1 = cell(1, nSubresultants);
arr_SingularValues = cell(1, nSubresultants);





% For each possible value of k, k = 0,...,min(m,n)
for i = 1 : 1 : nSubresultants
    
    % Set index i
    k = i + lowerLimit_k - 1;
    
    arr_Sk{i} = BuildSubresultant_3Polys(fx, gx, hx, k);
    
%     % If not the first subresultant, build by removing rows and columns.
%     if i > 1
%         
%         % update C_f and C_g by removing rows and columns
%         C_f1 = C_f1(1:m+n-k+1, 1:n-k+1);
%         C_f2 = C_f2(1:m+o-k+1, 1:o-k+1);
%         
%         C_g = C_g(1:m+n-k+1, 1:m-k+1);
%         C_h = C_h(1:m+o-k+1, 1:m-k+1);
%         
%         block = blkdiag(C_f1,C_f2);
%         col = [C_g ; C_h];
%         
%         arr_Sk{i} = [block col];
%         [~,R] = qr(arr_Sk{i});
%         
%         arr_R{i} = abs(R);
%         
%         %         % Perform QR Decomposition of Sk, by QR delete.
%         %         % Remove the last column
%         %         [Q,R] = qrdelete(Q,R,m+n+2-((2*k)-2),'col');
%         %
%         %         % Remove last column of C_{1}(f)
%         %         [Q,R] = qrdelete(Q,R,n+2-k,'col');
%         %
%         %         % Remove last row
%         %         [Q,R] = qrdelete(Q,R,m+n+2-k,'row');
%         
%     end
    
%     % Get number of rows in R1_{k}
%     [nRowsR1,~] = size(diag(arr_R{i}));
%     
%     % Obtain R1 the top square of the |R| matrix.
%     arr_R1{i} = arr_R{i}(1:nRowsR1, 1:nRowsR1);
    
    
    
    
end



global SETTINGS

%   * Singular Values
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Residuals

switch SETTINGS.METRIC
    
    case 'Minimum Singular Values'
        
        for i = 1:1:nSubresultants
            
            % Add to the vector of minimum Singular values from SVD of S_{k}.
            arr_SingularValues{i} = svd(arr_Sk{i});
            
            % Get the minimum Singular value from SVD of S_{k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
            
        end
        
        plotSingularValues(arr_SingularValues, limits_k, limits_t);
        %plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range);
        
        metric = log10(vMinimumSingularValues);
        
    case 'Residuals'
        
        for i = 1 : 1 : nSubresultants

            % Add to the vector of minimum distances
            vMinimumResidual(i) = GetMinDistance(arr_Sk);

        end
        
        metric = vMinimumResidual;
        
    case 'R1 Row Diagonals'
        
        
        for i = 1:1: nSubresultants
            
            vDiagsR1 = diag(arr_R1{i});
            
            % Get maximum and minimum row diagonals of R1
            vMaxDiagR1(i) = max(vDiagsR1);
            vMinDiagR1(i) = min(vDiagsR1);
        end
        metric = vMaxDiagR1./vMinDiagR1;
        
    case 'R1 Row Norms'
        
        for i = 1:1:nSubresultants
            
            % Get Norms of each row in the matrix R1
            vR1_RowNorms = sqrt(sum(arr_R1{i}.^2,2))./norm(arr_R1{i});
            
            % Get maximum and minimum row norms of rows of R1.
            vMaxRowNormR1(i) = max(vR1_RowNorms);
            vMinRowNormR1(i) = min(vR1_RowNorms);
            
        end
        
        metric = vMaxRowNormR1./vMinRowNormR1;
        
    otherwise
        error('err');
end

% If only one subresultant exists, use an alternative method.
if (upperLimit_k == lowerLimit_k ) % If only one Subresultant Exists
    
    if (lowerLimit_k == 1)
        % Use the singular values from the only subresultant S_{1} to determine
        % if S_{1} is full rank or rank deficient.
        
        t = GetGCDDegree_OneSubresultant(arr_SingularValues);
    else
        % Since can not be corpime, and only one subresultant exists
        
        t = lowerLimit_k;
        fprintf([mfilename ' : ' sprintf('One subresultant : t = %i \n',t)]);
    end
    
else
    % Get the type of problem.
    % Problem Type.
    % Singular      : All Subresultants S_{k} are Singular, and rank deficient
    % NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
    % Mixed         : Some Subresultants are Singular, others are Non-Singular.
    
    % PlotGraphs_3Polys()
    [t] = GetGCDDegree_MultipleSubresultants(metric, limits_k, limits_t, rank_range);
    
end





end






