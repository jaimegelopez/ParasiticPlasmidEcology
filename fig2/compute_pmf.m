function pmf = compute_pmf(beta,kappa)

%This function computes the WF pmf

pmf = zeros(length(beta),1);

pmf(1) = 1;

for ii = 2:length(pmf)
    pmf(ii) = compute_prod(beta,ii,kappa);
end

pmf = pmf./sum(pmf);

end

