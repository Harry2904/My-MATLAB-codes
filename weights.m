function w = weights(scheme, k)
    w = zeros(1, 3); 
    if scheme == 1
        w = [0, 1, 0];
    elseif scheme == 2
        if k == 2
            w = [0, 2, -1];
        else
            w = [0, 3/2, -1/2];
        end
    elseif scheme == 3
        if k == 2
            w = [1/3, 1, -1/3];
        else
            w = [3/8, 6/8, -1/8];
        end
    else
        error("Invalid scheme selected. Choose 1, 2, or 3.");
    end
end
