function [f_root_mult_arr, g_root_mult_arr, d_root_mult_arr, ...
    u_root_mult_arr, v_root_mult_arr] = GCD_Examples_Univariate_2Polys(ex_num)
%
%
%
% Inputs.
%
% ex_num : Example Number
%
% Outputs.
%
% f : Array of symbolic factors of f(x) and multiplicity of the factors.
%
% g : Array of symbolic factors of g(x) and multiplicity of the factors.
%
% d : Array of symbolic factors of d(x) and multiplicity of the factors.
%
% u : Array of symbolic factors of u(x) and multiplicity of the factors.
%
% v : Array of symbolic factors of v(x) and multiplicity of the factors.


syms x

switch ex_num
    
    
    
    case '1'
        
        d_root_mult_arr = [...
            (x + 1)   2
            (x + 2)   1
            ];
        u_root_mult_arr = [...
            (x + 3)   1
            ];
        v_root_mult_arr = [...
            (x - 2)   1
            ];
        
        
        
    case '2'
        % From The computation of the degree of an approximate greatest
        % common divisor of two Bernstein Polynomials - Bourne, Winkler
        % & Yi.
        
        u_root_mult_arr = [
            (x - 0.1)    1
            (x - 0.3)     2
            (x - 0.5)     2
            (x - 0.7)     3
            (x - 2.5)     3
            (x - 3.4)     3
            ];
        
        v_root_mult_arr = [
            (x - 0.85)    4
            (x - 0.9)     4
            (x - 1.1)     3
            ];
        
        d_root_mult_arr = [
            (x - 0.10)    3
            (x - 0.80)    2
            ];
        
   
        
    case '3'
        % From The computation of the degree of an approximate greatest
        % common divisor of two Bernstein Polynomials - Bourne, Winkler
        % & Yi.
        
        d_root_mult_arr = [...
            (x - 0.1)       2
            (x - 0.56)      4
            (x - 0.75)      3
            (x - 1.37)      3
            ];
        
        u_root_mult_arr = [
            (x - 0.1)       1
            (x - 0.82)      3
            (x + 0.27)      3
            (x - 1.46)      2
            ];
        
        v_root_mult_arr = [
            (x - 0.99)      4
            (x - 2.12)      1
            (x - 1.2)       3
            ];
        
    
    case '4'
        d_root_mult_arr = [...
            (x-0.5)         2
            ];
        u_root_mult_arr = [...
            (x + 1.234)     3
            ];
        v_root_mult_arr = [...
            (x-1.75292)     4
            ];
       
    case '5'
        
        d_root_mult_arr = [...
            (x - 1.2)         4
            (x + 4.7562)      2
            ];
        u_root_mult_arr = [...
            (x + 1.234)     3
            ];
        v_root_mult_arr = [...
            (x - 0.5)     3
            ];
      
    case '6'
        
        d_root_mult_arr = [...
            (x - 1.2)         14
            (x + 4.7562)      2
            ];
        u_root_mult_arr = [...
            (x - 1.5)     10
            ];
        v_root_mult_arr = [...
            (x - 0.5)     7
            (x + 1.9)     5
            ];
     
        
    case '7'
        
        d_root_mult_arr = [...
            (x - 1.2)         14
            (x + 4.7562)      2
            ];
        u_root_mult_arr = [...
            (x - 1.9)     10
            (x - 0.2)     4
            ];
        v_root_mult_arr = [...
            (x - 0.5)     7
            (x + 0.9)     15
            ];
      
        
    case '8'
        % From The computation of the degree of an approximate greatest
        % common divisor of two Bernstein Polynomials - Bourne, Winkler
        % & Yi.
        % ADAPTED
        
        d_root_mult_arr = [...
            (x - 0.1)       2
            (x - 0.56)      8
            (x - 0.75)      3
            (x - 1.37)      3
            ];
        
        u_root_mult_arr = [
            (x - 0.1)       1
            (x - 0.82)      3
            (x + 0.27)      3
            (x - 1.46)      2
            ];
        
        v_root_mult_arr = [
            (x - 0.99)      4
            (x - 2.12)      1
            (x - 1.2)       3
            ];
        
       
    case '9'
        % From The computation of the degree of an approximate greatest
        % common divisor of two Bernstein Polynomials - Bourne, Winkler
        % & Yi.
        % ADAPTED
        
        d_root_mult_arr = [...
            (x - 0.1)       2
            (x - 0.56)      8
            (x - 0.75)      10
            (x - 1.37)      3
            ];
        
        u_root_mult_arr = [
            (x - 0.1)       4
            (x - 0.82)      3
            (x + 0.27)      3
            (x - 1.46)      2
            ];
        
        v_root_mult_arr = [
            (x - 0.99)      4
            (x - 2.12)      1
            (x - 1.2)       3
            ];
        
       
        
        
        
    case '10'
        d_root_mult_arr = [...
            (x - 1.2435487954)       2
            (x - 5.56)      8
            (x - 0.75)      10
            (x - 1.37)      3
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        4
            (x - 0.82)              3
            (x + 2.27564657)        3
            (x - 1.46)              2
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)      4
            (x - 2.12)      1
            (x - 1.2222222)       3
            ];
       
        
    case '11'
        d_root_mult_arr = [...
            (x - 1.2435487954)      2
            (x - 5.56)              8
            (x - 0.7515678)         15
            (x - 1.37)              3
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        4
            (x - 0.82)              3
            (x + 2.27564657)        3
            (x - 1.46)              4
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        4
            (x - 2.12)              4
            (x - 1.2222222)         3
            ];
        
        
        
    case '12'
        d_root_mult_arr = [...
            (x - 1.2435487954)      1
            (x - 5.56)              2
            (x - 9.7515678)         1
            (x - 1.37)              1
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        1
            (x - 0.82)              1
            (x + 2.27564657)        2
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        1
            (x - 2.12)              1
            (x - 1.2222222)         1
            ];
       
    case '13'
        
        d_root_mult_arr = [...
            (x + 0.1)   2
            (x + 0.2)   1
            ];
        u_root_mult_arr = [...
            (x + 0.3)   1
            ];
        v_root_mult_arr = [...
            (x - 0.2)   1
            ];
        
       
        
    case '14'
        
        d_root_mult_arr = [...
            (x - 1.2435487954)      1
            (x - 0.5678910112)      7
            (x - 9.7515678)         5
            (x - 1.37)              3
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        1
            (x - 0.82)              1
            (x + 0.27564657)        10
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        9
            (x - 2.12)              1
            (x - 1.2222222)         1
            ];
        
       
        
    case '15'
        d_root_mult_arr = [...
            (x - 1.2435487954)      1
            (x + 5.103579)          2
            (x - 0.76549843)        4
            (x - 0.37987984)        1
            ];
        
        u_root_mult_arr = [
            (x - 0.5465444984)      1
            (x - 0.7165792)         1
            (x + 2.27564657)        2
            (x - 1.600548798)       1
            ];
        
        v_root_mult_arr = [
            (x - 0.9916546598)      1
            (x - 0.54987)           1
            (x - 1.229876852)       1
            ];
      
        
        
        
    case '16'
        d_root_mult_arr = [...
            (x - 0.7515678)         7
            (x - 0.37)              1
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        1
            (x - 0.82)             30
            (x + 2.27564657)        4
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        1
            (x - 0.12)              1
            (x - 1.2222222)         1
            ];
      
        
    case '17'
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              1
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        1
            (x - 0.82)             30
            (x + 2.27564657)        2
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        1
            (x - 0.12)              1
            (x - 1.2222222)         1
            ];
      
    case '18'
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              1
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        1
            (x - 0.82)              30
            (x + 2.27564657)        20
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        1
            (x - 0.12)              1
            (x - 1.2222222)         1
            ];
       
        
    case '19' % m >> n >> t
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              5
            ];
        
        u_root_mult_arr = [
            (x - 0.10122344)        1
            (x - 0.82)              30
            (x + 2.27564657)        20
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        5
            (x - 0.12)              5
            (x - 1.2222222)         1
            ];
       
        
        
    case '20' 
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              5
            ];
        
        u_root_mult_arr = [
            (x - 1.46)              1
            (x - 0.82)              10
            (x - 0.10122344)        1
            (x + 2.27564657)        20
            
            ];
        
        v_root_mult_arr = [
            (x - 1.2222222)         1
            (x - 0.99102445)        10
            (x - 0.12)              20
            ];
       
        
    case '21' 
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              5
            ];
        
        u_root_mult_arr = [
            (x - 0.82)              10
            (x + 2.27564657)        20
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        10
            (x - 1.2222222)         1
            (x - 0.1)               20
            ];
       
        
    case '22' 
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              5
            ];
        
        u_root_mult_arr = [
            (x - 0.82)              10
            (x + 2.27564657)        20
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        10
            (x - 1.2222222)         1
            (x - 1.1)               20
            ];
        
        
        
         case '22a' 
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              5
            ];
        
        u_root_mult_arr = [
            (x - 0.82)              1
            (x + 2.27564657)        2
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        1
            (x - 1.2222222)         1
            (x - 1.1)               2
            ];
        
        
           case '22b' 
        d_root_mult_arr = [...
            (x - 0.7515678)         2
            (x - 0.37)              5
            ];
        
        u_root_mult_arr = [
            (x - 0.82)              10
            (x + 2.27564657)        10
            (x - 1.46)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        10
            (x - 1.2222222)         1
            (x - 1.1)               10
            ];
        
    case '23' % case 15 edited
        d_root_mult_arr = [...
            (x - 1.2435487954)      1
            (x + 1.56)              3
            (x - 0.2176387)         3
            (x - 0.157981)          1
            ];
        
        u_root_mult_arr = [
            (x - 0.65487987)        1
            (x - 0.82657984)        5
            (x + 0.27564657)        2
            (x - 5.26)              1
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        1
            (x - 0.12)              8
            (x - 0.4579879)         1
            ];
        
        
       case '24' % case 15 edited
        d_root_mult_arr = [...
            (x - 1.2435487954)      5
            (x + 1.56)              7
            (x - 0.217612343)         3
            (x - 0.157981)          9
            ];
        
        u_root_mult_arr = [
            (x - 0.85487987)        1
            (x - 5.8265747784)      7
            (x + 0.27564657)        5
            (x - 0.26)              2
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        5
            (x + 0.12)              4
            (x - 0.4579879)         9
            
            ];  
        case '25' % case 15 edited
        d_root_mult_arr = [...
            (x - 1.2435487954)      1
            (x + 1.56)              7
            (x - 0.217612343)       3
            (x - 0.157981)          9
            ];
        
        u_root_mult_arr = [
            (x - 0.85487987)        1
            (x - 0.8265747784)        5
            (x + 0.27564657)        2
            (x - 0.26)              2
            ];
        
        v_root_mult_arr = [
            (x - 0.99102445)        5
            (x + 0.12)              4
            (x - 4.4579879)         1
            ];  
        
        
        
    case 'Random 1'
        
        d_root_mult_arr = BuildRandomPoly(10, -1, 1, 1, 5);
        u_root_mult_arr = BuildRandomPoly(4, -1, 1, 1, 5);
        v_root_mult_arr = BuildRandomPoly(7, -1, 1, 1, 5);
        
        
       
        
    case 'Random 2'
        
        d_root_mult_arr = BuildRandomPoly(3, -1, 1, 1, 5);
        u_root_mult_arr = BuildRandomPoly(7, -2, 1, 1, 7);
        v_root_mult_arr = BuildRandomPoly(4, -1, 1, 1, 5);
        
        
    
        
        
        
        
%         case '7'
%         
%         d_root_mult_arr = [...
%             (x - 1.2)         14
%             (x + 4.7562)      2
%             ];
%         u_root_mult_arr = [...
%             (x - 1.9)     10
%             (x - 0.2)     4
%             ];
%         v_root_mult_arr = [...
%             (x - 0.5)     7
%             (x + 0.9)     15
%             ];
        
   
        case 'test7a'
        
        d_root_mult_arr = [...
            (x - 1.2)         1
            (x + 4.7562)      2
            ];
        u_root_mult_arr = [...
            (x - 1.9)     1
            (x - 0.2)     4
            ];
        v_root_mult_arr = [...
            (x - 0.5)     3
            (x + 0.9)     1
            ];
        
          case 'test7b'
        
        d_root_mult_arr = [...
            (x - 1.2)         2
            (x + 4.7562)      4
            ];
        u_root_mult_arr = [...
            (x - 1.9)     1
            (x - 0.2)     4
            ];
        v_root_mult_arr = [...
            (x - 0.5)     3
            (x + 0.9)     1
            ];
        
            case 'test7c'
        
        d_root_mult_arr = [...
            (x - 1.2)         2
            (x + 4.7562)      4
            ];
        u_root_mult_arr = [...
            (x - 1.9)     2
            (x - 0.2)     8
            ];
        v_root_mult_arr = [...
            (x - 0.5)     6
            (x + 0.9)     2
            ];
        
          case 'test7d'
        
        d_root_mult_arr = [...
            (x - 1.2)         4
            (x + 2.7562)      8
            ];
        u_root_mult_arr = [...
            (x - 0.9)     4
            (x - 0.2)     15
            ];
        v_root_mult_arr = [...
            (x - 1.5)     6
            (x + 0.9)     2
            ];

        
    case {'1a', '1b', '1c'} 
        [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num);    
              
    case {'2a', '2b', '2c'} 
        [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num);
    
    case {'3a', '3b', '3c'} 
        [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num);
        
    case {'4a', '4b', '4c'} 
        [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num);
    
    case {'10a', '10b', '10c'} 
        [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num);
        
        
    case {'11a', '11b', '11c'}
        [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num);
        
    case {'200a', '200b', '200c'}
        [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num);
        
        
        
    otherwise
        
        pattern = 'Random:m=(\d+).n=(\d+).t=(\d+).low=(-?\d+).high=(-?\d+).seed=(\d+)';
        
        if ~isempty(regexp(ex_num,pattern,'start'))
            
            str = ex_num;
            
            expression_m = regexp(str,'m=(\d+)','tokens');
            m_str = expression_m{1};
            
            expression_n = regexp(str,'n=(\d+)','tokens');
            n_str = expression_n{1};
            
            expression_t = regexp(str,'t=(\d+)','tokens');
            t_str = expression_t{1};
            
            expression_low = regexp(str, 'low=(-?\d+)', 'tokens');
            low_str = expression_low{1};
            
            expression_high = regexp(str, 'high=(-?\d+)','tokens');
            high_str = expression_high{1};
            
            expression_seed = regexp(str, 'seed=(\d+)','tokens');
            seed_str = expression_seed{1};
            
        end
        
        
        
        
        m = str2double(m_str);
        n = str2double(n_str);
        t = str2double(t_str);
        interval_low = str2double(low_str);
        interval_high = str2double(high_str);
        seed = str2double(seed_str);
        
        rng(seed)
        
        
        m_t = m - t;
        n_t = n - t;
        
        
        % Get a vector of multiplicities
        desired_sum = t;
        nRoots_dx = randi([1, t ], 1, 1); % Choose the number of roots
        vMult_dx = diff([0,sort(randperm(desired_sum - 1, nRoots_dx - 1)), desired_sum]);
        
        
        desired_sum = m_t; % Choose the desired sum
        nRoots_ux = randi([1, m_t ], 1, 1); % Choose the number of integers
        vMult_ux = diff([0,sort(randperm(desired_sum - 1, nRoots_ux - 1)), desired_sum]);
        
        desired_sum = n_t; % Choose the desired sum
        nRoots_vx = randi([1, n_t ], 1, 1); % Choose the number of integers
        vMult_vx = diff([0,sort(randperm(desired_sum - 1, nRoots_vx - 1)), desired_sum]);
        
        
        nRoots = length(vMult_dx) + length(vMult_ux) + length(vMult_vx);
        
        % Generate set of distinct roots
        vRandom = floor(rand(nRoots, 1)*1e10)/1e10;
        vRoots = interval_low + vRandom * (interval_high - interval_low);
        
        d_root_mult_arr = [ x - vRoots(1:nRoots_dx) vMult_dx'];
        
        u_root_mult_arr = [x - vRoots(nRoots_dx + 1: nRoots_dx + nRoots_ux) vMult_ux'];
        
        v_root_mult_arr = [x - vRoots(nRoots_dx + nRoots_ux + 1 : end) vMult_vx'];
        
        f_root_mult_arr = [u_root_mult_arr; d_root_mult_arr];
        g_root_mult_arr = [v_root_mult_arr; d_root_mult_arr];
        
        
        
        
        
        
end

f_root_mult_arr = [u_root_mult_arr ; d_root_mult_arr];
g_root_mult_arr = [v_root_mult_arr ; d_root_mult_arr];


end

function [d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ThreePolyProblem(ex_num)

example_number = ex_num(1 : end - 1);
example_variant = ex_num(end);

[~, ~, ~, ...
    d_rm, u_rm, v_rm, w_rm] = GCD_Examples_Univariate_3Polys(example_number);

switch example_variant
    
    case 'a'
        d_root_mult_arr = d_rm;
        u_root_mult_arr = u_rm;
        v_root_mult_arr = v_rm;
        
    case 'b'
        d_root_mult_arr = d_rm;
        u_root_mult_arr = u_rm;
        v_root_mult_arr = w_rm;
        
    case 'c'
        d_root_mult_arr = d_rm;
        u_root_mult_arr = v_rm;
        v_root_mult_arr = w_rm;
end

end


function [root_mult_arr] = BuildRandomPoly(n_factors, root_low, root_high, mult_low, mult_high)
%
% % Inputs
%
% n_factors : (Int) Number of factors
%
% root_low : (Float)
%
% root_high : (Float)
%
% mult_low : (Int)
%
% mult_high : (Int)

syms x

root_mult_arr = sym(zeros(n_factors, 2));

for i = 1 : 1 : n_factors
    
    a = root_low;
    b = root_high;
    root = a + (b-a).*rand(1,1);
    
    factor = x - root;
    
    
    
    mult = randi([mult_low mult_high],1,1);
    
    root_mult_arr(i,:) = [factor ; mult];
    
end



end