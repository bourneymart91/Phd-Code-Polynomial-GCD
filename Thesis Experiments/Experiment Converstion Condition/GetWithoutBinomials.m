function [f] = GetWithoutBinomials(f_bi)

m = GetDegree(f_bi);
bi_m = GetBinomials(m);

f = f_bi ./ bi_m;

end