%% Set-up
clear;clc

%This script creates a version of the phase diagrams with finite loss rates

%Parameters
par.Gamma = 1;
par.delta1 = 0.1;
par.alpha = 1;
par.delta2 = 0.3;
par.Np = 0.6;
par.n = 1000;
par.rl = 5e-2;
par.alphaval = 0.8;
Delta = logspace(-0.7,-0.01,par.n);
transfer_rates = logspace(-4,1,par.n);

%Make grids
[Delta_mat,transfer_mat] = meshgrid(Delta,transfer_rates);


%% Make transformation diagram

%Generate symbolic functions
[jacobian_cell,sol,~] = symbolic_transf_loss();

sol1 = [sol.c(1),sol.rho(1),sol.rhop(1),sol.p(1)];
sol2 = [sol.c(2),sol.rho(2),sol.rhop(2),sol.p(2)];
sol3 = [sol.c(3),sol.rho(3),sol.rhop(3),sol.p(3)];

sol_fcn = cell(3,1);
sol_fcn{1} = matlabFunction(sol1);
sol_fcn{2} = matlabFunction(sol2);
sol_fcn{3} = matlabFunction(sol3);

jacobian_fcn = cell(3,1);
jacobian_fcn{1} = matlabFunction(jacobian_cell{1});
jacobian_fcn{2} = matlabFunction(jacobian_cell{2});
jacobian_fcn{3} = matlabFunction(jacobian_cell{3});

sol_stability{1} = nan(par.n,par.n);
sol_stability{2} = nan(par.n,par.n);
sol_stability{3} = nan(par.n,par.n);

sol_phys = sol_stability;

%Numerically generate stability from the Jacobian
for i = 1:par.n
    for j = 1:par.n
        for k = 1:3
            par.gammat = transfer_mat(i,j);
            par.Delta = Delta_mat(i,j);
            
            if k == 1
                numeric_eig = eig(jacobian_fcn{k}(par.Delta,par.Gamma,par.Np,par.alpha,par.delta1,par.delta2,par.gammat,par.rl));
                numeric_sol = sol_fcn{k}(par.Gamma,par.alpha,par.delta1);
            else
                numeric_eig = eig(jacobian_fcn{k}(par.Delta,par.Gamma,par.Np,par.alpha,par.delta1,par.delta2,par.gammat,par.rl));
                numeric_sol = sol_fcn{k}(par.Delta,par.Gamma,par.Np,par.alpha,par.delta1,par.delta2,par.gammat,par.rl);
            end
            
            sol_stability{k}(i,j) = (max(real((numeric_eig))) < 0) & (min(numeric_sol) >= 0) & isreal(numeric_sol);
            sol_phys{k}(i,j) = (min(numeric_sol) >= 0) & isreal(numeric_sol);
            
            if k == 2
                sol2_rho(i,j) = numeric_sol(2);
            elseif k == 1
                sol1_rho(i,j) = numeric_sol(2);
            else
                sol3_rho(i,j) = numeric_sol(2);
            end
            
        end
    end
end

%Generate the full diagram
total_trans_grid = zeros(par.n,par.n);
total_trans_grid(sol_stability{1} == 1) = 1;
total_trans_grid(sol_stability{2} == 1) = 2;
total_trans_grid((sol_stability{1} == 1) & (sol_stability{2} == 1)) = 3;
total_trans_grid((sol_stability{1} ~= 1) & (sol_stability{2} ~= 1)) = 4;

third_stability_trans = sum(sum(sol_stability{3}));
agreement_trans = sum(~(total_trans_grid(:) == 4) & ~sol_stability{3}(:))/sum(~sol_stability{3}(:));

%% Make conjugation diagram

%Generate symbolic functions
[jacobian_cell,sol,sol_labels] = symbolic_conj_loss();

sol1 = [sol.c(1),sol.rho(1),sol.rhop(1)];
sol2 = [sol.c(2),sol.rho(2),sol.rhop(2)];
sol3 = [sol.c(3),sol.rho(3),sol.rhop(3)];

sol_fcn = cell(3,1);
sol_fcn{1} = matlabFunction(sol1);
sol_fcn{2} = matlabFunction(sol2);
sol_fcn{3} = matlabFunction(sol3);

jacobian_fcn = cell(3,1);
jacobian_fcn{1} = matlabFunction(jacobian_cell{1});
jacobian_fcn{2} = matlabFunction(jacobian_cell{2});
jacobian_fcn{3} = matlabFunction(jacobian_cell{3});

sol_stability{1} = nan(par.n,par.n);
sol_stability{2} = nan(par.n,par.n);
sol_stability{3} = nan(par.n,par.n);

sol_phys = sol_stability;

%Loops through and numerically compute Jacobian eigenvalues
for i = 1:par.n
    for j = 1:par.n
        for k = 1:3
            par.gammac = transfer_mat(i,j);
            par.Delta = Delta_mat(i,j);
                if k == 1
                    numeric_eig = eig(jacobian_fcn{k}(par.Delta,par.Gamma,par.alpha,par.delta1,par.gammac,par.rl));
                    numeric_sol = sol_fcn{k}(par.Gamma,par.alpha,par.delta1);
                else
                    numeric_eig = eig(jacobian_fcn{k}(par.Delta,par.Gamma,par.alpha,par.delta1,par.gammac,par.rl));
                    numeric_sol = sol_fcn{k}(par.Delta,par.Gamma,par.alpha,par.delta1,par.gammac,par.rl);
                end

            sol_stability{k}(i,j) = (max(real((numeric_eig))) < 0) & (min(numeric_sol) >= 0) & isreal(numeric_sol);
            sol_phys{k}(i,j) = (min(numeric_sol) >= 0)& isreal(numeric_sol);
            
            if k == 2
                sol2_rho(i,j) = numeric_sol(2);
            elseif k == 1
                sol1_rho(i,j) = numeric_sol(2);
            else
                sol3_rho(i,j) = numeric_sol(2);
            end
            
        end
    end
end


%Generate phase diagram
total_free = zeros(par.n,par.n);
total_free(sol_stability{1} == 1) = sol1_rho(sol_stability{1} == 1);
total_free(sol_stability{2} == 1) = sol2_rho(sol_stability{2} == 1);
total_free(sol_stability{3} == 1) = sol3_rho(sol_stability{3} == 1);

third_stability_conj = sum(sum(sol_stability{3}));
multi_stable = sum((sol_stability{1}(:) + sol_stability{2}(:) + sol_stability{3}(:)) > 1);

save('phase_with_loss_data.mat');

%% Plot transformation diagram

load('phase_with_loss_data.mat');

newfigure(3.42/2, (3/3)*3.42/2);

gap = [0.1,0.15];
marg_h = [0.08,0.04];
marg_w = [0.12,0.34];
Fig1ax = tight_subplot(2,1,gap,marg_h,marg_w);

%Order is No plasmid, plasmid only, bistable, coexistence
mymap = [1 0.85 0; 0.9 0.425 0; 1 1 1; 0.8 0 0];
colormap(Fig1ax(2),mymap);

axes(Fig1ax(2));

%Create checkboard grid of states 1 and 2 to impose over the bistability
%(state 3)
bistable_pattern = double(checkerboard(50,10,10) > 0.5);
bistable_pattern(bistable_pattern == 0) = 2;
total_trans_grid(total_trans_grid == 3) = bistable_pattern(total_trans_grid == 3);

%Plot
hold on
im2 = imagesc(total_trans_grid,[1,4]);
im2.AlphaData = par.alphaval;

xlabel({'$\log$(Plasmid cost, $\Delta$)'},'Interpreter','latex')
ylabel({'$\log$(Trans. rate, $\gamma_\textrm{t}$)'},'Interpreter','latex')

xlim([1,par.n]);

ylim([1,par.n]);

yticks([])
xticks([])
box on
set(gca,'FontSize',4)

text(-0.1,1.05,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)

%% Plot conjugation diagram

axes(Fig1ax(1));
colormap(parula)
imagesc(flipud(total_free)./(par.Gamma/par.delta1))
h = colorbar;
caxis([0,1]);
h.Position = h.Position + [0.35,0,0,0];
ylabel(h,'Fraction of plasmid-free cells','Interpreter','latex');
h.Ticks = [0,0.5,1];
h.TickLabelInterpreter = 'latex';
xlabel({'$\log$(Plasmid cost, $\Delta$)'},'Interpreter','latex')
ylabel({'$\log$(Conj. rate, $\gamma_\textrm{c}$)'},'Interpreter','latex')

xlim([1,par.n]);

ylim([1,par.n]);

yticks([])
xticks([])
box on
set(gca,'FontSize',4)

text(-0.1,1.05,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)


%% Make legend

axes(Fig1ax(2));
ax = gca;    
ax.Clipping = 'off';   
level = 0.8*par.n;
rec_width = 0.05*par.n;
make_rec = @(coord,color)rectangle('Position',...
    [coord(1)-rec_width/2,coord(2)-rec_width/2,rec_width,rec_width],...
    'FaceColor',color,'EdgeColor',color);
make_half_rec = @(coord,color)rectangle('Position',...
    [coord(1)-rec_width/4,coord(2)-rec_width/2,rec_width/2,rec_width],...
    'FaceColor',color,'EdgeColor',color);

point_offset = 0.05*par.n;
point1 = [1.1*par.n,level];
point2 = point1 + [0,-150];
make_rec(point2,mymap(1,:))
text(point2(1) + point_offset,point2(2),'No plasmid',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point3 = point2 + [0,-150];
make_rec(point3,mymap(2,:))
text(point3(1) + point_offset,point3(2),'Coexistence',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point4 = point3 + [0,-150];
make_half_rec(point4 - [rec_width/4,0],mymap(2,:))
make_half_rec(point4 + [rec_width/4,0],mymap(1,:))
text(point4(1) + point_offset,point4(2),'Bistability',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

print(gcf, '-dpng','supp_figure_phase_with_loss.png','-r600','-painters');

close all