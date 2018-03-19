function [PolyA_rm_array, PolyB_rm_array, PolyC_rm_array, ...
    PolyD_rm_array, PolyP_rm_array, PolyQ_rm_array, PolyR_rm_array] ...
    = GCD_Examples_Univariate_3Polys(ex_num)
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% % Outputs.
%
% PolyA_rm_array : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial a(x)
%
% PolyB_rm_array  : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial b(x)
%
% PolyC_rm_array : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial c(x)
%
% PolyD_rm_array : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial d(x)
%
% PolyP_rm_array : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial p(x)
%
% PolyQ_rm_array : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial q(x)
%
% PolyR_rm_array : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial r(x)

syms x


% 
[PolyD_rm_array, PolyP_rm_array, PolyQ_rm_array, PolyR_rm_array] = ...
    GetPolys(ex_num);


% Get the factors and corresponding multiplicities of the polynomials a(x),
% b(x) and c(x) which are the polynomials f(x), g(x) and h(x) where the 
% ordering has not yet been defined.

PolyA_rm_array = [PolyP_rm_array;  PolyD_rm_array];
PolyB_rm_array = [PolyQ_rm_array;  PolyD_rm_array];
PolyC_rm_array = [PolyR_rm_array;  PolyD_rm_array];






end




function [PolyD_rm_arr, PolyP_rm_arr, PolyQ_rm_arr, PolyR_rm_arr] = GetPolys(ex_num)
%
% % Inputs
%
%
% ex_num : (String) Example number
%
%
%
% % Outputs
%
%
% PolyD_rm_arr : (Matrix) Symbolic factors and corresponding multiplicities of the
% polynomial d(x). Where d(x) is the GCD of f(x), g(x) and h(x).
%
% PolyP_rm_arr : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial P(x)
%
% PolyQ_rm_arr : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial Q(x)
%
% PolyR_rm_arr : (Matrix) Matrix containing the symbolic factors and 
% corresponding multiplicities of the polynomial R(x)

syms x


% Initialise matrices to store common factors (and corresponding
% muliplicities) of the pairwise polynomials.
Common_PolyP_PolyQ = ...
    [
    ];

Common_PolyQ_PolyR = ...
    [
    ];

Common_PolyP_PolyR = ...
    [
    ];







switch ex_num
    case '1'
        
        PolyD_rm_arr = ...
            [...
            (x - 0.1)   2
            (x - 0.2)   1
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 0.3)   1
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 0.2)   1
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.512)   2
            ];
        
        
        
        
        
    case '2'
        
        
        PolyP_rm_arr = [
            (x - 0.1)    1
            (x - 0.3)     2
            (x - 0.5)     2
            (x - 0.7)     3
            (x - 2.5)     3
            (x - 3.4)     3
            ];
        
        PolyQ_rm_arr = [
            (x - 0.85)    4
            (x - 0.9)     4
            (x - 1.1)     3
            ];
        
        PolyR_rm_arr = ...
            [...
            (x + 1.75)    4
            ];
        
        PolyD_rm_arr = [
            (x - 0.10)    3
            (x - 0.80)    5
            ];
        
    case '2Variant'
        
        
        PolyP_rm_arr = [
            (x - 0.1)    1
            (x - 0.3)     2
            (x - 0.5)     2
            (x - 0.7)     3
            (x - 0.5)     3
            (x - 0.4)     3
            ];
        
        PolyQ_rm_arr = [
            (x - 0.85)    4
            (x - 0.9)     4
            (x - 0.1)     3
            ];
        
        PolyR_rm_arr = ...
            [...
            (x + 0.75)    4
            ];
        
        PolyD_rm_arr = [
            (x - 0.10)    3
            (x - 0.80)    5
            ];
        
        
        
        
    case '3'
        
        
        PolyD_rm_arr = [...
            (x - 0.1)       2
            (x - 0.56)      4
            (x - 0.75)      5
            (x - 1.37)      3
            ];
        
        PolyP_rm_arr = [
            (x - 0.82)      3
            (x + 0.27)      4
            (x - 1.46)      2
            ];
        
        PolyQ_rm_arr = [
            (x - 0.99)      4
            (x - 2.12)      4
            
            ];
        
        PolyR_rm_arr = ...
            [...
            (x + 1.75)      2
            (x - 5.72)      8
            ];
        
        Common_PolyP_PolyQ = ...
            [
            (x - 1.2)       3
            ];
        
    case '3Variant'
        
        
        PolyD_rm_arr = [...
            (x - 0.1)       2
            (x - 0.56)      4
            (x - 0.75)      5
            (x - 1.37)      3
            ];
        
        PolyP_rm_arr = [
            (x - 0.82)      3
            (x + 0.27)      4
            (x - 1.46)      2
            ];
        
        PolyQ_rm_arr = [
            (x - 0.99)      4
            (x - 2.12)      4
            
            ];
        
        PolyR_rm_arr = ...
            [...
            (x + 1.75)      2
            (x - 5.72)      8
            ];
        
        Common_PolyP_PolyQ = ...
            [
            (x - 1.2)       6
            ];
        
    case '4'
        PolyD_rm_arr = ...
            [...
            (x - 0.5)         2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x + 1.234)     3
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75292)     4
            ];
        
        PolyR_rm_arr = ...
            [
            (x + 1.75)  2
            ];
        
        
        
    case '5'
        
        PolyD_rm_arr = ...
            [...
            (x - 1.2)         4
            (x + 4.7562)      2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x + 1.234)     3
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 0.5)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x + 1.75)  2
            ];
        
    case '6'
        PolyD_rm_arr = ...
            [...
            (x - 0.5654654561)      5
            (x - 0.21657894)        5
            ];
        
        PolyP_rm_arr = ...
            [...
            (x + 1.234)             3
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75292)           4
            (x - 0.99851354877)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x + 1.75)              2
            (x + 0.778912324654)    4
            ];
        
        
        
    case '7'
        PolyD_rm_arr = ...
            [...
            (x - 0.5654654561)      5
            (x - 0.21657894)        10
            (x + 1.2468796514)      3
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 0.234)             3
            (x - 0.01564897)        2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75292)           4
            (x - 0.99851354877)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 1.75)              2
            (x - 0.778912324654)    4
            ];
        
        
    case '8'
        PolyD_rm_arr = ...
            [...
            (x - 0.5654654561)      5
            (x - 0.21657894)        1
            (x + 0.2468796514)      3
            ];
        
        PolyP_rm_arr = ...
            [...
            (x + 0.7879734)             1
            (x - 0.01564897)        2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 0.75292)           20
            (x - 0.99851354877)     7
            ];
        
        PolyR_rm_arr = ...
            [
            (x + 1.75)              2
            (x + 0.778912324654)    4
            ];
        
    case '8Variant'
        PolyD_rm_arr = ...
            [...
            (x - 0.56)      5
            (x - 0.21657894)        1
            (x + 0.246)      3
            ];
        
        PolyP_rm_arr = ...
            [...
            (x + 0.78)             1
            (x - 0.01)        2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 0.75)           5
            (x + 0.75)           5
            (x - 0.9)     7
            ];
        
        PolyR_rm_arr = ...
            [
            (x + 1.75)              2
            (x + 0.778912324654)    4
            ];
        
        
    case '9'
        PolyD_rm_arr = ...
            [...
            (x - 0.5654654561)      2
            (x - 0.21657894)        1
            ];
        
        PolyP_rm_arr = ...
            [...
            (x + 0.7879734)         6
            (x - 0.01564897)        6
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 0.75292)           1
            (x - 0.99851354877)     1
            ];
        
        PolyR_rm_arr = ...
            [
            (x + 1.75)              2
            (x + 0.778912324654)    2
            ];
        
    case '10'
        PolyD_rm_arr = ...
            [...
            (x - 0.5654654561)      2
            (x - 0.21657894)        1
            (x + 1.654987654)       4
            (x - 1.2657984335)      4
            ];
        
        PolyP_rm_arr = ...
            [...
            (x + 0.7879734)         3
            (x - 0.41564897)        6
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 0.75292)           1
            (x - 0.99851354877)     4
            (x + 0.1654988136)      4
            ];
        
        PolyR_rm_arr = ...
            [
            (x + 1.75)              2
            (x - 0.564987986958)    3
            (x + 0.778912324654)    2
            ];
        
    case '11'
        PolyD_rm_arr = ...
            [...
            (x - 9.2657984335)      2
            (x - 1.2657984335)      4
            (x - 0.21657894)        1
            (x - 0.0654654561)      2
            (x + 1.654987654)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 0.41564897)        6
            (x + 0.7879734)         9
            (x + 1.932654987)       1
            (x + 2.3549879)         2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75292)           1
            (x - 0.99851354877)     3
            (x + 0.1654988136)      4
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.564987986958)    3
            (x + 0.778912324654)    2
            (x + 1.75)              2
            ];
    case '12'
        
        
        PolyD_rm_arr = [...
            (x - 0.1)       2
            (x - 0.56)      4
            (x - 0.75)      5
            (x - 1.37)      3
            ];
        
        PolyP_rm_arr = [
            (x - 1.2)       5
            (x - 0.82)      5
            (x + 0.27)      10
            (x - 1.46)      10
            ];
        
        PolyQ_rm_arr = [
            (x - 0.99)      10
            (x - 2.12)      10
            (x - 1.2)       10
            ];
        
        PolyR_rm_arr = ...
            [...
            (x + 1.75)      2
            (x - 5.72)      4
            ];
        
        
    case '14'
        
        PolyD_rm_arr = [...
            (x - 0.1)       10
            (x - 0.56)      4
            (x - 0.75)      5
            (x - 1.37)      3
            ];
        
        PolyP_rm_arr = [
            (x - 1.2)       1
            (x - 0.82)      3
            (x + 0.27)      4
            (x - 1.46)      20
            ];
        
        PolyQ_rm_arr = [
            (x - 0.99)      4
            (x - 2.12)      4
            (x - 1.2)       13
            ];
        
        PolyR_rm_arr = ...
            [...
            (x + 1.75)      2
            (x - 5.72)      12
            ];
        
    case '15'
        
        PolyD_rm_arr = ...
            [...
            (x - 0.1)   10
            (x - 0.2)   1
            ];
        PolyP_rm_arr = ...
            [...
            (x - 0.3)   20
            ];
        PolyQ_rm_arr = ...
            [...
            (x - 0.2)   5
            ];
        PolyR_rm_arr = ...
            [
            (x - 0.512)   7
            ];
        
        
    case '16' % Modified version of 3  where roots are of similar magnitude, so coefficients are of similar magnitude
        
        
        
        PolyD_rm_arr = [...
            (x - 0.1)       2
            (x - 0.56)      4
            (x - 0.75)      5
            (x - 1.37)      3
            ];
        
        PolyP_rm_arr = [
            (x - 1.2)       1
            (x - 0.82)      3
            (x + 0.27)      4
            (x - 1.46)      2
            ];
        
        PolyQ_rm_arr = [
            (x - 0.99)      4
            (x - 0.12)      4
            (x + 0.2)       3
            ];
        
        PolyR_rm_arr = ...
            [...
            (x + 0.75)      2
            (x - 0.72)      8
            ];
        
        
        
    case '100'
        PolyD_rm_arr = ...
            [...
            (x - 0.2)      4
            (x - 0.06)      2
            (x + 1.65)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 0.41)        4
            (x - 0.35)         2
            (x + 0.78)         2
            (x + 1.93)       1
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75)           1
            (x - 0.99)     3
            (x + 0.16)      4
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.48)    3
            (x + 0.71)    4
            (x + 1.75)              2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            (x - 0.11 )     4
            (x - 0.34)        1
            (x - 0.56)      2
            ];
        
        Common_PolyQ_PolyR = ...
            [
            (x - 0.22) 3
            ];
        
        Common_PolyP_PolyR = ...
            [
            ];
        
    case '101'
        PolyD_rm_arr = ...
            [...
            (x - 0.2)      1
            (x - 0.06)      2
            (x + 1.65)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 0.41)        4
            (x - 0.35)         2
            (x + 0.78)         2
            (x + 1.93)       1
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75)           1
            (x - 0.99)     3
            (x + 0.16)      4
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.48)    3
            (x + 0.71)    4
            (x + 1.75)              2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            ( x - 0.65)    2
            ( x - 0.23)    4
            ];
        
        Common_PolyQ_PolyR = ...
            [
            ];
        
        Common_PolyP_PolyR = ...
            [
            ];
        
        
    case '102'
        PolyD_rm_arr = ...
            [...
            (x - 0.2654987)      1
            (x - 0.06654)      1
            (x + 1.65654)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 1.41)        4
            (x - 0.35)         2
            (x + 0.78)         2
            (x + 1.93)       1
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75)           1
            (x - 0.99)     3
            (x + 0.16)      4
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.48)    3
            (x - 0.71)    4
            (x + 0.75)    2
            (x - 1.2)     2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            ( x - 0.65)    2
            ( x - 0.23)    6
            ];
        
        Common_PolyQ_PolyR = ...
            [
            ( x + 1.82)    2
            ( x - 0.79)    3
            ];
        
        Common_PolyP_PolyR = ...
            [
            ];
        
        
    case '103'
        PolyD_rm_arr = ...
            [...
            (x - 0.2654987)      4
            (x - 1.65654)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 1.41)        4
            (x - 0.35)         2
            (x + 1.93)       2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 1.75)           1
            (x - 0.99)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.71)    4
            (x + 0.75)    2
            (x - 1.2)     2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            ( x - 0.65)    2
            ( x - 7.23)    3
            ];
        
        Common_PolyQ_PolyR = ...
            [
            (x - 0.1234)      1
            (x - 0.56789)     1
            ];
        
        Common_PolyP_PolyR = ...
            [
            
            ];
        
    case '104'
        PolyD_rm_arr = ...
            [...
            (x - 0.2654987)      4
            (x - 1.65654)       3
            (x + 12.45 )        2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 1.41)        4
            (x - 0.35)         2
            (x + 1.93)       2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 0.75)           1
            (x - 7.99)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.71)    4
            (x + 0.75)    2
            (x - 1.2)     2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            ( x + 0.65)    2
            ( x - 7.23)    3
            ];
        
        Common_PolyQ_PolyR = ...
            [
            (x - 0.1234)      1
            (x - 0.56789)     3
            ];
        
        Common_PolyP_PolyR = ...
            [
            (x - 0.92357)      2
            (x + 1.56789)     4
            ];
        
        
        
        
        
        
        
        
        
        %Examples for use in power basis
        
        
    case '200'
        PolyD_rm_arr = ...
            [...
            (x - 1.26)      4
            (x - 4.65)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 1)         4
            (x - 1.5)       2
            (x + 3)         2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 2)           1
            (x - 4.99)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.71)    4
            (x + 1.75)    2
            (x - 3.2)     2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            (x - 1.1)    3
            ];
        
        Common_PolyQ_PolyR = ...
            [
            ];
        
        Common_PolyP_PolyR = ...
            [
            ];
        
    case '201'
        PolyD_rm_arr = ...
            [...
            (x - 1.57)      4
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 1)         4
            (x - 1.5)       2
            (x + 3)         2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 2.9)           1
            (x - 4.99)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 0.71)    4
            (x + 1.75)    2
            (x - 1.2)     2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            (x - 0.5)  4
            ];
        
        Common_PolyQ_PolyR = ...
            [
            ];
        
        Common_PolyP_PolyR = ...
            [
            ];
        
    case '202'
        PolyD_rm_arr = ...
            [...
            (x - 1.2)      4
            (x - 4.6)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 1)         4
            (x - 1.5)       2
            (x + 3)         2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 2)           1
            (x - 0.9)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 3.71)    4
            (x + 0.75)    2
            (x - 1.2)     2
            (x - 1.8)     2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            (x - 0.5)  5
            ];
        
        Common_PolyQ_PolyR = ...
            [
            ];
        
        Common_PolyP_PolyR = ...
            [
            ];
        
    case '203'
        PolyD_rm_arr = ...
            [...
            (x - 0.1)       2
            ];
        
        PolyP_rm_arr = ...
            [...
            (x - 1.9)         4
            (x + 3)         2
            ];
        
        PolyQ_rm_arr = ...
            [...
            (x - 2.1)           1
            (x + 0.9)     3
            ];
        
        PolyR_rm_arr = ...
            [
            (x - 3.71)    4
            (x + 0.75)    2
            (x - 1.2)     2
            (x - 1.8)     2
            ];
        
        
        Common_PolyP_PolyQ = ...
            [
            (x - 0.5)  5
            ];
        
        Common_PolyQ_PolyR = ...
            [
            (x - 1.7)   3
            ];
        
        Common_PolyP_PolyR = ...
            [
            ];
        
        
end

% Add the pairwise common factors to matrices of factors of P(x), Q(x) 
% and R(x)
PolyP_rm_arr = [PolyP_rm_arr ; Common_PolyP_PolyQ; Common_PolyP_PolyR ];
PolyQ_rm_arr = [PolyQ_rm_arr ; Common_PolyP_PolyQ; Common_PolyQ_PolyR ];
PolyR_rm_arr = [PolyR_rm_arr ; Common_PolyP_PolyR; Common_PolyQ_PolyR ];




end
