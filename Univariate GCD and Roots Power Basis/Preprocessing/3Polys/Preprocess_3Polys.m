function [GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta] = Preprocess_3Polys(fx,gx,hx)

global SETTINGS

switch SETTINGS.N_EQUATIONS_SYLVESTER_MATRIX

    case '2'
    
        [GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta] = Preprocess_3Polys_2Eqns(fx, gx, hx);


    case '3'
        
        [GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta] = Preprocess_3Polys_3Eqns(fx, gx, hx);
   
        
end


end