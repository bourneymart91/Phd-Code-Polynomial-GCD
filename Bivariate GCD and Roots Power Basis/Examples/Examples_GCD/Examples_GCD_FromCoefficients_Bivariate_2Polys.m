function [fxy, gxy,uxy,vxy,dxy,...
    m,m1,m2,...
    n,n1,n2,...
    t,t1,t2] = Examples_GCD_FromCoefficients_Bivariate_2Polys(ex_num)


syms x y

addpath(genpath(('../Examples')))
[f_roots_mult_arr,g_roots_mult_arr,d_roots_mult_arr,...
    u_roots_mult_arr,v_roots_mult_arr] = GCD_Examples_Bivariate_2Polys(ex_num);



[fxy,m,m1,m2] = GetCoefficientsFromSymbolicRoots(f_roots_mult_arr);
[gxy,n,n1,n2] = GetCoefficientsFromSymbolicRoots(g_roots_mult_arr);
[dxy,t,t1,t2] = GetCoefficientsFromSymbolicRoots(d_roots_mult_arr);
[uxy,~,~,~] = GetCoefficientsFromSymbolicRoots(u_roots_mult_arr);
[vxy,~,~,~] = GetCoefficientsFromSymbolicRoots(v_roots_mult_arr);

symbolic_f = GetSymbolicPolyFromSymbolicRoots(f_roots_mult_arr);
symbolic_g = GetSymbolicPolyFromSymbolicRoots(g_roots_mult_arr);
symbolic_d = GetSymbolicPolyFromSymbolicRoots(d_roots_mult_arr);

display(symbolic_f);
display(symbolic_g);
display(symbolic_d);



end