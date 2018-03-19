function [f_noisy]=Noise(f,el,eu)

% Add noise to the coefficients of polynomial f, given in matrix form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% f :   Matrix of coefficients of polynomial f(x,y)

% el :  Signal to noise low limit

% eu :  Signal to noise upper limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the degree of input polynomial f
m1 = size(f,1) - 1;
m2 = size(f,2) - 1; 

switch nargin
    case 2 % Only one noise is specified, set upper = lower
               
        rp = (2*rand(m1+1,m2+1))-ones(m1+1,m2+1);
        s = rp*el;
        
        noise_matrix = f.*s;
        f_noisy = f + noise_matrix;
        
        
    case 3 % Specified upper and lower bound of noise
        
        y = (2*rand(m1+1,m2+1))-ones(m1+1,m2+1);
        s = eu *ones(m1+1,m2+1) -  y.*(eu-el);
        noise_matrix = f'.*s;
        f_noisy = f + noise_matrix';
end