function w_n = weights_nonuniform(scheme, Ds_D, Ds_U, Ds_UU)
    w_n = zeros(1, 3);
    
    if scheme == 1
        w_n(1) = 0;
        w_n(2) = 1;
        w_n(3) = 0;
        
    elseif scheme == 2
        w_n(1) = Ds_U * (2 * Ds_U + Ds_UU);
        w_n(1) = w_n(1) / ((Ds_D + Ds_U) * (Ds_D + 2 * Ds_U + Ds_UU));
        
        w_n(2) = Ds_D * (2 * Ds_U + Ds_UU);
        w_n(2) = w_n(2) / ((Ds_D + Ds_U) * (Ds_U + Ds_UU));
        
        w_n(3) = -Ds_D * Ds_U;
        w_n(3) = w_n(3) / ((Ds_U + Ds_UU) * (Ds_D + 2 * Ds_U + Ds_UU));
    end
end