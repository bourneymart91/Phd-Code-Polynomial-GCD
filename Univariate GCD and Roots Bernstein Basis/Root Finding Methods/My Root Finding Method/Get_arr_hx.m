function arr_hx = Get_arr_hx(arr_fx)

global SETTINGS

switch SETTINGS.GET_ARR_HX_METHOD
    
    case 'From Deconvolution'
    
        fprintf([mfilename ' : ' sprintf('Deconvolving f_{i}(x) : %s \n',SETTINGS.DECONVOLUTION_METHOD_HX_FX)])
        arr_hx = Deconvolve_Set(arr_fx,SETTINGS.DECONVOLUTION_METHOD_HX_FX);
    
    case 'From ux'
        
        arr_hx = arr_hx;
        
    otherwise
        error('err');
end