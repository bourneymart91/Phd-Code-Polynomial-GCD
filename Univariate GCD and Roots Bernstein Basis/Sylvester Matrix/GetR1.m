function R1 = GetR1(Sk)
% Get the square upper triangular matrix R1 from the QR Decomposition of
% the kth subresultant matrix Sk
%
% Inputs.
%
% Sk : k-th Subresultant matrix S_{k}(f,g)

% Get QR Decomposition of the sylvester matrix S_{k}(f(\theta,w),g(\theta,w))
[~,R] = qr(Sk);

% Take absolute values of R_{k}
R = abs(R);

% Get number of rows in R1_{k}
[nRowsR1,~] = size(diag(R));

% Obtain R1 the top square of the R matrix.
R1 = R(1:nRowsR1,1:nRowsR1);


end