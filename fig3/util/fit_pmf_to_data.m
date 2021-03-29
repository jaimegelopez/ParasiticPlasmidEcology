function [fitted_pmf, fitted_par] = fit_pmf_to_data(data_pmf,exp_vec)

%This fits a pmf to the data with a prespecified exponent form

Delta0 = 0.01;
q0 = 0.005;
trunc = length(data_pmf);

x0 = [Delta0,q0];

fitted_par = lsqnonlin(@(x) theory_data_err(x,exp_vec,data_pmf,trunc), x0);

fitted_pmf = get_theory_pmf(fitted_par(1),fitted_par(2),exp_vec);

end

