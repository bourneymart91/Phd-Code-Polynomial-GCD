function [wx,wy,wxy] = o_roots_mymethod_xy(wx,wy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform a series of GCD calculations on the w_{x,i}s

LineBreakLarge()

% Get number of w_{x,i}
[~,nEntries_wx] = size(wx);

wxy = {};

for i = 1:1:nEntries_wx
    
    % Get number of columns in w{i}(x)
    [nRows,nCols] = size(wx{i});
    
    % If nCols > 1, then factor is bivariate
    if nCols > 1
        
        % Check that w_{i}(y) has more than one row
        if nRows > 1
            
            fxy = wx{i};
            gxy = wy{i};
            alpha = 1;
            th1 = 1;
            th2 = 1;
            
            [~,nCols] = size(fxy);
            [nRows,~] = size(gxy);
            t1 = nRows -1;
            t2 = nCols -1;
                      
            
            t = t1;
            % Get Quotients
            [uxy_matrix_clc,vxy_matrix_clc] = GetCofactors(fxy,gxy,t);
            
            % Get the GCD dxy
            dxy_matrix_clc = GetGCDCoefficients(fxy,gxy,uxy_matrix_clc,vxy_matrix_clc,alpha, th1, th2);

            
            %Overwrite wx and wy with new values
            wxy{i} = dxy_matrix_clc;
            wx{i} = uxy_matrix_clc;
            wy_new = Deconvolve_Bivariate(wy{i},dxy_matrix_clc);
            wy{i} = vxy_matrix_clc;
            
            fprintf([mfilename ' : ' sprintf('Roots of degree %i',i)])
            display(wxy{i})
            display(wx{i})
            display(wy{i})
        else
            error('err')
        end
        
        
    end
end


