clear;clc

%This script looks at convergence to the binomial distribution

%% Declare variables

syms rho0 rho1 c gammac Gamma alpha1 f delta1 rl

%% Solve one-plasmid

one1 = alpha1*c*rho0 - gammac*rho0*rho1 - delta1*rho0 + alpha1*rl*(1-f)*c*rho1;
one2 = (1 - f)*alpha1*c*rho1 + gammac*rho0*rho1 - delta1*rho1 - alpha1*rl*(1-f)*c*rho1;
one3 = Gamma - alpha1*c*rho0 - (1-f)*alpha1*c*rho1; 

one_vars = [rho0,rho1,c];

one_eqns = [one1 == 0, one2 == 0, one3 == 0];

one_sol = solve(one_eqns,one_vars);

%% Find which is the coexistence equilibria

par.delta1 = 1;
par.f = 1e-2;
par.rl = 0.05;
par.alpha1 = 1;
par.Gamma = 1;
gammac_crit = (par.delta1*(par.f) + par.delta1*par.rl*(1-par.f))*(par.delta1/par.Gamma);
par.gammac = 3*gammac_crit;
par.n_plasmid = 3;

rho0_one_fcn =  matlabFunction(one_sol.rho0);
rho1_one_fcn = matlabFunction(one_sol.rho1);

rho0_one = rho0_one_fcn(par.Gamma,par.alpha1,par.delta1,par.f,par.gammac,par.rl);
rho1_one = rho1_one_fcn(par.Gamma,par.alpha1,par.delta1,par.f,par.gammac,par.rl);

%Find the coexistence equilibria
rho0_ind = (rho0_one > 0) & (rho0_one < par.Gamma/par.delta1);
rho1_ind = (rho1_one > 0) & (rho1_one < par.Gamma/par.delta1);
coex_ind = find(rho0_ind & rho1_ind);

coex_rho0 = one_sol.rho0(coex_ind);
coex_rho1 = one_sol.rho1(coex_ind);

%% Build binomial solution for three plasmids

one_p = coex_rho1/(coex_rho0 + coex_rho1);

bino_rho0 = (Gamma/delta1)*1*(one_p^0)*(1-one_p)^3;
bino_rho1 = (Gamma/delta1)*3*(one_p^1)*(1-one_p)^2;
bino_rho2 = (Gamma/delta1)*3*(one_p^2)*(1-one_p)^1;
bino_rho3 = (Gamma/delta1)*1*(one_p^3)*(1-one_p)^0;

%% Compare bino to three plasmid numbers

f_vec = logspace(-5,-1,100);

rho0_bino_fcn = matlabFunction(bino_rho0);
rho1_bino_fcn = matlabFunction(bino_rho1);
rho2_bino_fcn = matlabFunction(bino_rho2);
rho3_bino_fcn = matlabFunction(bino_rho3);

real_sol_mat = zeros(length(f_vec),4);
bino_sol_mat = zeros(length(f_vec),4);
bino_sol_error = zeros(length(f_vec),1);
bino_RHS_error = zeros(length(f_vec),1);
ref_norm = zeros(length(f_vec),1);

par.additive_fitness = 0;
par.error_threshold = 1e-7;
par.dt = 0.1*min([1/par.alpha1,1/par.gammac]);


for i = 1:length(f_vec)

    par.f = f_vec(i);
    gammac_crit = (par.delta1*(par.f) + par.delta1*par.rl/par.alpha1);
    par.gammac = 3*gammac_crit;
    
    rho0_numeric_bino = rho0_bino_fcn(par.Gamma,par.alpha1,par.delta1,par.f,par.gammac,par.rl);
    rho1_numeric_bino = rho1_bino_fcn(par.Gamma,par.alpha1,par.delta1,par.f,par.gammac,par.rl);
    rho2_numeric_bino = rho2_bino_fcn(par.Gamma,par.alpha1,par.delta1,par.f,par.gammac,par.rl);
    rho3_numeric_bino = rho3_bino_fcn(par.Gamma,par.alpha1,par.delta1,par.f,par.gammac,par.rl);
    
    bino_sol_mat(i,:) = [rho0_numeric_bino,rho1_numeric_bino,rho2_numeric_bino,rho3_numeric_bino];
    
    [counts] = simulate_condensed_plasmid(par);
    
    real_sol_mat(i,:) = counts';
    
    bino_sol_error(i) = norm(real_sol_mat(i,:) - bino_sol_mat(i,:));
    
end

%% Plot error vs cost

newfigure(0.75*3.42/2, 0.75*(2/3)*3.42/2);

loglog(f_vec,bino_sol_error,'k-','LineWidth',0.7)
xlabel('Plasmid cost, $\Delta$','Interpreter','latex')
ylabel({'Deviation from','binomial solution'},'Interpreter','latex')
xlim([1e-5,1e-1])
xticks([1e-5,1e-4,1e-3,1e-2,1e-1])
yticks([1e-6,1e-4,1e-2,1e0])
ylim([1e-6,1e0])
set(gca,'FontSize',4);
set(gca,'TickLabelInterpreter','latex')
axest = gca;
axest.XAxis.MinorTick = 'off';
axest.YAxis.MinorTick = 'off';
box off
print(gcf, '-dpng','supp_figure_binomial_convergence.png','-r600','-painters');