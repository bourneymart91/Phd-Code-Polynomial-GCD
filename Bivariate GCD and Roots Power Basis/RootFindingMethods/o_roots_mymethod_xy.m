function [wx,wy,wxy] = o_roots_mymethod_xy(wx,wy,vDegt_wx,vDegt_wy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform a series of GCD calculations on the w_{x,i}s

LineBreakLarge()

% Get number of w_{x,i}
[nPolys_arr_wx] = size(wx,1);

wxy = {};

for i = 1:1:nPolys_arr_wx
    
    % Get number of columns in w{i}(x)
    [~,nCols] = size(wx{i});
    
    [nRows,~] = size(wy{i});
    % If nCols > 1, then factor is bivariate
    if nCols > 1
        
        % Check that w_{i}(y) has more than one row
        if nRows > 1
            
            fxy_matrix_n = wx{i};
            gxy_matrix_n = wy{i};
            alpha = 1;
            th1 = 1;
            th2 = 1;
            
            [~,nCols] = size(fxy_matrix_n);
            [nRows,~] = size(gxy_matrix_n);
            t1 = nRows -1;
            t2 = nCols -1;
            
            % % Get Quotients
            
            % Preprocess f and g
            fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);
            gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);
            
            
            [uww_matrix_clc,vww_matrix_clc] = GetQuotients_Relative(fww_matrix,alpha.*gww_matrix,t1,t2);
            
            % Get the GCD d(x,y)
            dxy_matrix_clc = GetGCDCoefficients(fww_matrix,alpha.*gww_matrix,uww_matrix_clc,vww_matrix_clc,m,n,t);

            
            %Overwrite wx and wy with new values
            wxy{i} = dxy_matrix_clc;
            wx{i} = uww_matrix_clc;
            wy_new = Deconvolve_Bivariate_Single_Respective(wy{i},dxy_matrix_clc);
            wy{i} = vww_matrix_clc;
            
            fprintf([mfilename ' : ' sprintf('Roots of degree %i',i)])
            display(wxy{i})
            display(wx{i})
            display(wy{i})
        else
            error('err')
        end
        
        
    end
end


