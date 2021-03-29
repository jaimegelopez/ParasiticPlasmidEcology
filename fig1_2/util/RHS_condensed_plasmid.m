function dy = RHS_condensed_plasmid(y,gamma,rl,cN,cK,Dv,Fv,alpha,Gamma,cond_combos,num_variants)

%Computes the RHS for the condensed plasmid model

rhovec = y(1:end-1);
c = y(end);

%Get vectors for matrix operations
single_rho = rhovec./num_variants;
combo_rho = (rhovec(cond_combos(:,1)+1)./num_variants(cond_combos(:,1)+1)).*...
    (rhovec(cond_combos(:,2)+1)./num_variants(cond_combos(:,2)+1)); 

drho = alpha.*Fv.*rhovec.*c + gamma*cN*combo_rho + alpha.*c.*rl.*cK*single_rho - Dv.*rhovec;

dc = Gamma - sum(alpha.*Fv.*rhovec.*c); 

dy = [drho; dc];

end

