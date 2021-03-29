function [fitted_pmf, fitted_par] = fit_var_epis_pmf(data_pmf,base_exp_vec)

%This function look for an optimal fit by varying all parameters. Delta and q 
%searched by lsqnonlin. Epistasis exponent varied by grid search due to 
%lsqnonlin not handling the parameter well

Delta0 = 0.01;
q0 = 0.005;
trunc = length(data_pmf);
x0 = [Delta0,q0];

exp_range = linspace(0,2,100);

for i = 1:length(exp_range)
    [fitted_par_cell{i},err_vec(i)] = lsqnonlin(@(x) var_epis_err(x,exp_range(i),base_exp_vec,data_pmf,trunc), x0);
end

[~,min_ind] = min(err_vec);
fitted_par = fitted_par_cell{min_ind};
fitted_par(3) = exp_range(min_ind);

if fitted_par(3) == 0
    base_exp_vec(2:end) = 1;
    fitted_pmf = get_theory_pmf(fitted_par(1),fitted_par(2),base_exp_vec);
else
    fitted_pmf = get_theory_pmf(fitted_par(1),fitted_par(2),base_exp_vec.^exp_range(min_ind));
end

end

