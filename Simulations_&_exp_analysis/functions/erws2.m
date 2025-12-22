function result = erws2(ampa_current, gaba_current, params, v_0, rec_step)
    if nargin < 4
        v_0 = 0;
    end
    if nargin < 5
        rec_step = 2.0;
    end

    tau_AMPA = params(1) * power(v_0, -params(2)) + params(3);
    tau_GABA = params(4) * power(v_0, -params(5)) + params(6);
    alpha = params(7) * power(v_0, -params(8)) + params(9);

    ampa_delayed = circshift(ampa_current, round(tau_AMPA / rec_step));
    gaba_delayed = circshift(gaba_current, round(tau_GABA / rec_step));

    result = ampa_delayed - alpha * gaba_delayed;
    
end


