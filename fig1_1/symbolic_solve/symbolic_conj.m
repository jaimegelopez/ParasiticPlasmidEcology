function [jacobian_cell,sol,sol_labels] = symbolic_conj()

syms rho rhop c gammac Gamma alpha1 Delta delta1

base1 = alpha1*c*rho - gammac*rho*rhop - delta1*rho;
base2 = (1-Delta)*alpha1*c*rhop + gammac*rho*rhop - delta1*rhop;
base3 = Gamma - alpha1*c*rho - (1-Delta)*alpha1*c*rhop; 

vars = [rhop,c,rho];

eqns = [base1 == 0, ...
    base2 == 0,...
    base3 == 0];

sol = solve(eqns,vars);

jacobian_cell = {};
var_jacobian = jacobian([base2, base3, base1],vars);

sol_labels = cell(3,1);

for i = 1:length(sol.rho)
    jacobian_cell{i} = subs(var_jacobian,[rho,rhop,c],[sol.rho(i),sol.rhop(i),sol.c(i)]);
    if ~isequal(sol.rho(i),sym(0)) && ~isequal(sol.rhop(i),sym(0))
        sol_labels{i} = 'coexistence';
    elseif isequal(sol.rho(i),sym(0)) && ~isequal(sol.rhop(i),sym(0))
        sol_labels{i} = 'plasmid-only';
    else
        sol_labels{i} = 'plasmid-free';
    end
end



end