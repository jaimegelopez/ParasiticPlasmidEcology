function pmf = get_theory_pmf(Delta,q,exp_vec)

beta = (1-Delta).^(exp_vec);
kappa = q/(1-q);

pmf = zeros(length(beta),1);

pmf(1) = 1;

for ii = 2:length(pmf)
    pmf(ii) = compute_prod(beta,ii,kappa);
end

pmf = pmf./sum(pmf);

end

