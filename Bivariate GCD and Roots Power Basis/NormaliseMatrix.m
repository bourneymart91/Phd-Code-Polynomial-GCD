function fxy_n = NormaliseMatrix(fxy)

if fxy(1,1) ~= 0
    fxy_n = fxy./fxy(1,1);
else
    fxy_n = fxy;
end


end