function [T_nplus, T_nminus] = Temp_f_d(k, w_np, w_nm, T1, T2, T3, T4)
  
    if k < 1 || k > size(w_np,1)
        error('Index k is out of bounds for w_np and w_nm');
    end

    T_nplus = w_np(k,1)*T3 + (w_np(k,2) - 1)*T2 + w_np(k,3)*T1;
    T_nminus = w_nm(k,1)*T2 + (w_nm(k,2) - 1)*T3 + w_nm(k,3)*T4;
end

