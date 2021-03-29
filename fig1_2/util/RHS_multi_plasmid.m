function dy = RHS_multi_plasmid(y,gamma,rl,N,K,Dv,Fv,alpha,Gamma,type_combo_index)

%Computes the RHS for the explicit multi-plasmid model

%Get populations
sv = y(1:end-1);

%Get all possible conjugation combinations
vv = sv(type_combo_index(:,1)).*sv(type_combo_index(:,2));

%Get nutrient
c = y(end);

dsv = alpha.*Fv.*sv.*c + gamma*N*vv + alpha.*c.*rl.*K*sv - Dv.*sv;

dc = Gamma - sum(alpha.*Fv.*sv.*c);

dy = [dsv; dc];

end

