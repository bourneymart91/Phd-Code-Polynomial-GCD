function [lambda] = GetMean_3Polys_3Eqns(fx, n_k, o_k)
% Calculate Geometric means of input polynomials f(x,y) and g(x,y)


global SETTINGS


switch SETTINGS.MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        lambda = GetGeometricMeanMatlabMethod_3Polys_3Eqns(fx, n_k, o_k);
        
    case 'Geometric Mean My Method'
        lambda = GetGeometricMean_3Polys_3Eqns(fx, n_k);
        
        
    case 'Arithmetic Mean'
        
        lambda = GetArithmeticMean_3Polys_3Eqns(fx, n_k, m_k);
        
        
    case 'None'
        
        lambda = 1;
             
    otherwise
        error('err: MEAN_METHOD must be either *Geometric Mean Matlab Method or *Geometric Man My Method or *None')
end
end