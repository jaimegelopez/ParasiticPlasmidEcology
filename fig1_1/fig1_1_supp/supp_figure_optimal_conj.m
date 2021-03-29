clear;clc

%This script computes the optimal plasmid fitness costs and loss rates

%Set up parameters for sweep
par.max_pcn = 100;
par.min_pcn = 1.1;
par.n = 100;
par.Gamma = 1;
par.delta1 = 1;
par.alpha = 1;
par.loss_ymax = 1.1;

%Get vector of plasmid copy costs
Deltap_vec = logspace(-5,-2,par.n);
opt_pcns = zeros(1,par.n);
true_opt_pcns = opt_pcns;
loss_fit_ratio = opt_pcns;

%Loops through plasmid copy costs
for i = 1:par.n
    
par.Deltap = Deltap_vec(i);

cost_fun = @(np) par.delta1*(np.*par.Deltap + 2.^(1-np).*(1-np.*par.Deltap));

%Find optimal pcn
opt_pcns(i) = fminsearch(@(x) cost_fun(x),5);

min_loss = par.delta1*(2.^(1-opt_pcns(i)).*(1-opt_pcns(i).*par.Deltap));
min_fit = par.delta1*(opt_pcns(i).*par.Deltap);

%Find loss to cost ratio
loss_fit_ratio(i) = min_loss/min_fit;

true_opt_pcns(i) = (par.Deltap*(-lambertw(0,exp(1)*2^(1/par.Deltap - 1))) + par.Deltap + log(2))/(par.Deltap*log(2));

end

%Plot
figure
semilogx(Deltap_vec,opt_pcns,'k-')

%Plot supp figure
newfigure(0.75*3.42/2, 0.75*(2/3)*3.42/2);
semilogx(Deltap_vec,loss_fit_ratio,'k-','LineWidth',0.7)
xlabel('Plasmid copy cost, $\Delta_\textrm{p}$','Interpreter','latex')
ylabel('Ratio of optimal $p_\ell$ to $\Delta$','Interpreter','latex')
set(gca,'FontSize',4);
set(gca,'TickLabelInterpreter','latex')
axest = gca;
axest.XAxis.MinorTick = 'off';
axest.YAxis.MinorTick = 'off';
box off

print(gcf, '-dpng','supp_figure_optimal_conj.png','-r600','-painters');
