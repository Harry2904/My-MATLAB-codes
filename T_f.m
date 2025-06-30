function [Tf_plus, Tf_minus] = T_f(w, T1, T2, T3, T4)
    Tf_plus = w(1) * T3 + w(2) * T2 + w(3) * T1;
    Tf_minus = w(1) * T2 + w(2) * T3 + w(3) * T4;
end
