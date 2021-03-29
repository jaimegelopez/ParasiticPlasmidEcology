function err = theory_data_err(x,exp_vec,data_pmf,trunc)

%Compute error between theory and data

Delta = x(1);
q = x(2);

theory_pmf = get_theory_pmf(Delta,q,exp_vec);

err = norm(log(theory_pmf(1:trunc)) - log(data_pmf)');

end

