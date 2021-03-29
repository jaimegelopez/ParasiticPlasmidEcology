function [err,theory_pmf] = var_epis_err(x,epis_exp,base_exp_vec,data_pmf,trunc)

%Computes the error for a pmf with arbitrary epistasis exponent

Delta = x(1);
q = x(2);

if epis_exp == 0
    base_exp_vec(2:end) = 1;
    theory_pmf = get_theory_pmf(Delta,q,base_exp_vec);
else
    theory_pmf = get_theory_pmf(Delta,q,base_exp_vec.^epis_exp);
end

err = norm(log(theory_pmf(1:trunc)) - log(data_pmf(1:end)'));

end

