function [jacobian_cell,sol,sol_labels] = symbolic_transf()

syms rho rhop c gammat Gamma alpha1 Delta delta1 Np p delta2

base1 = alpha1*c*rho - gammat*rho*p - delta1*rho;
base2 = (1-Delta)*alpha1*c*rhop + gammat*rho*p - delta1*rhop;
base3 = Gamma - alpha1*c*rho - (1-Delta)*alpha1*c*rhop; 
base4 = Np*delta1*rhop - gammat*rho*p - delta2*p;

base_expr = [base1,base3,base2,base4];
vars = [rho,c,rhop,p];

eqns = [base1 == 0, base2 == 0, base3 == 0, base4 == 0];

sol = solve(eqns,vars);

jacobian_cell = {};
var_jacobian = jacobian(base_expr,vars);
sol_labels = cell(3,1);
for i = 1:length(sol.rho)
    jacobian_cell{i} = subs(var_jacobian,[rho,rhop,c,p],[sol.rho(i),sol.rhop(i),sol.c(i),sol.p(i)]);
    if ~isequal(sol.rho(i),sym(0)) && ~isequal(sol.rhop(i),sym(0))
        sol_labels{i} = 'coexistence';
    elseif isequal(sol.rho(i),sym(0)) && ~isequal(sol.rhop(i),sym(0))
        sol_labels{i} = 'plasmid-only';
    else
        sol_labels{i} = 'plasmid-free';
    end
end


end