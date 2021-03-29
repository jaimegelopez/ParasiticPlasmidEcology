function [product] = compute_prod(beta,ii,kappa)

%This function computes products in the pmf sum

product = 1;
for jj = 2:ii
    product = kappa.*product.*(beta(jj-1)./(1 - beta(jj)));
end


end
