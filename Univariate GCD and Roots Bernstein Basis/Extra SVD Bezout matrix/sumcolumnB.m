
function [] = sumcolumnB(B,fignum)

% This functions calculates the sum of the absolute values of the
% elements in each column of the Bezout matrix B=B(f,g), where
% f and g are polynomials.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum the absolute value of the elements in each column.
absB=abs(B);
colsum=sum(absB);
logcolsum=log10(colsum);

i=1:size(B,1);
figure(fignum)
plot(i,logcolsum,'b-o','MarkerSize',6,'MarkerFaceColor','b','LineWidth',1)
xlabel('column')
ylabel('log_{10} column sum')