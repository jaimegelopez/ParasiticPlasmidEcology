function [agreement_conj,agreement_transf] = fig1_1_phase_diagram(Fig1ax,par,loc,labels)

%% Set-up

rho_star = par.Gamma/par.delta1; 
Delta = logspace(-0.7,-0.01,par.n);
transfer_rates = logspace(-4,1,par.n);
%Make grids
[Delta_mat,transfer_mat] = meshgrid(Delta,transfer_rates);

%Order is No plasmid, plasmid only, bistable, coexistence
%mymap = [0 0 1; 1 0 0; 1 1 1; 1 0 1];
mymap = [1 0.85 0; 0.8 0 0; 1 1 1; 0.9 0.425 0];
colormap(mymap);


%% Generate the conjugation phase diagram

axes(Fig1ax(loc(1)));

%Define the stability functions
conj_plasmid_free_stable = @(Delta) Delta.*par.delta1./rho_star;
conj_plasmid_only_stable = @(Delta) ((1./(1-Delta))-1)*par.delta1./rho_star;

%Get the phase boundaries
conj_plasmid_free_gamma = conj_plasmid_free_stable(Delta);
conj_plasmid_only_gamma = conj_plasmid_only_stable(Delta);

%Get the mesh grid identities
conj_plasmid_free_grid = transfer_mat < conj_plasmid_free_stable(Delta_mat);
conj_plasmid_only_grid = transfer_mat > conj_plasmid_only_stable(Delta_mat);

%Build total grid
total_conj_grid = zeros(par.n,par.n);
total_conj_grid(conj_plasmid_free_grid == 1 & conj_plasmid_only_grid == 0) = 1;
total_conj_grid(conj_plasmid_free_grid == 0 & conj_plasmid_only_grid == 1) = 2;
total_conj_grid(conj_plasmid_free_grid == 1 & conj_plasmid_only_grid == 1) = 3;
total_conj_grid(conj_plasmid_free_grid == 0 & conj_plasmid_only_grid == 0) = 4;

%Look at coexistence stability
[jacobian_cell,sol,sol_labels] = symbolic_conj();
coex_index = strcmp(sol_labels,'coexistence');
coex_sol = [sol.c(coex_index),sol.rho(coex_index),sol.rhop(coex_index)];
sol_fcn = matlabFunction(coex_sol);
coex_jacobian = jacobian_cell{coex_index};
coex_eig = eig(coex_jacobian);
eig_fcn = matlabFunction(coex_eig);
coex_stability = zeros(size(transfer_mat));

for i = 1:par.n
    for j = 1:par.n
        par.gammac = transfer_mat(i,j);
        par.Delta = Delta_mat(i,j);
        numeric_eig = eig_fcn(par.Delta,par.Gamma,par.alpha,par.delta1,par.gammac);
        numeric_sol = sol_fcn(par.Delta,par.Gamma,par.alpha,par.delta1,par.gammac);
        coex_stability(i,j) = (max(real((numeric_eig))) < 0) & (min(numeric_sol) > 0) & isreal(numeric_sol);
    end
end

predicted_coex = total_conj_grid == 4;
agreement_conj = sum(sum(predicted_coex == 1 & coex_stability == 1))/sum(sum(coex_stability==1));

hold on
im1 = imagesc(total_conj_grid,[1,4]);
im1.AlphaData = par.alphaval;
plot(norm2im(Delta,Delta,par.n),norm2im(conj_plasmid_free_gamma,transfer_rates,par.n),'k-')
plot(norm2im(Delta,Delta,par.n),norm2im(conj_plasmid_only_gamma,transfer_rates,par.n),'k-')

xlabel({'$\log$(Plasmid cost, $\Delta$)'},'Interpreter','latex')
ylabel({'$\log$(Conjug.','rate, $\gamma_\textrm{c}$)'},'Interpreter','latex')

xlim([1,par.n]);
ylim([1,par.n]);

yticks([])
xticks([])
box on
set(gca,'FontSize',4)
text(-0.25,1.12,labels{1},'Interpreter','latex','Units','normalized','FontSize',4)

%Make legend
ax = gca;    
ax.Clipping = 'off';   
level = -0.7*par.n;
rec_width = 0.05*par.n;
make_rec = @(coord,color)rectangle('Position',...
    [coord(1)-rec_width/2,coord(2)-rec_width/2,rec_width,rec_width],...
    'FaceColor',color,'EdgeColor',color);
make_half_rec = @(coord,color)rectangle('Position',...
    [coord(1)-rec_width/4,coord(2)-rec_width/2,rec_width/2,rec_width],...
    'FaceColor',color,'EdgeColor',color);

point_offset = 0.05*par.n;
point1 = [-0.15*par.n,level];
make_rec(point1,mymap(2,:))
text(point1(1) + point_offset,point1(2),'Plasmid only',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point2 = point1 + [0.75*par.n,0];
make_rec(point2,mymap(1,:))
text(point2(1) + point_offset,point2(2),'No plasmid',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point3 = point2 + [0.69*par.n,0];
make_rec(point3,mymap(4,:))
text(point3(1) + point_offset,point3(2),'Coexistence',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point4 = point3 + [0.7*par.n,0];
make_half_rec(point4 - [rec_width/4,0],mymap(2,:))
make_half_rec(point4 + [rec_width/4,0],mymap(1,:))
text(point4(1) + point_offset,point4(2),'Bistability',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

%% Generate transformation phase diagram

axes(Fig1ax(loc(2)))

trans_plasmid_free_stable = @(Delta) par.delta2.*Delta./(par.Np -Delta)*(1./rho_star);
trans_plasmid_only_stable = @(Delta) par.delta1.*((1./(1 - Delta)) - 1).*(par.delta2./(par.Gamma.*par.Np));

trans_plasmid_free_gamma = trans_plasmid_free_stable(Delta);
trans_plasmid_only_gamma = trans_plasmid_only_stable(Delta);

%Get the mesh grid identities
trans_plasmid_free_grid = transfer_mat < trans_plasmid_free_stable(Delta_mat) | par.Np < Delta;
trans_plasmid_only_grid = transfer_mat > trans_plasmid_only_stable(Delta_mat);

%Build total grid
total_trans_grid = zeros(par.n,par.n);
total_trans_grid(trans_plasmid_free_grid == 1 & trans_plasmid_only_grid == 0) = 1;
total_trans_grid(trans_plasmid_free_grid == 0 & trans_plasmid_only_grid == 1) = 2;
total_trans_grid(trans_plasmid_free_grid == 1 & trans_plasmid_only_grid == 1) = 3;
total_trans_grid(trans_plasmid_free_grid == 0 & trans_plasmid_only_grid == 0) = 4;

%Look at coexistence stability
[jacobian_cell,sol,sol_labels] = symbolic_transf();
coex_ind = strcmp(sol_labels,'coexistence');
coex_sol = [sol.c(coex_ind),sol.rho(coex_ind),sol.rhop(coex_ind),sol.p(coex_ind)];
sol_fcn = matlabFunction(coex_sol);
coex_jacobian = jacobian_cell{coex_ind};
jacobian_fcn = matlabFunction(coex_jacobian);
coex_stability = zeros(size(transfer_mat));
for i = 1:par.n
    for j = 1:par.n
        par.gammat = transfer_mat(i,j);
        par.Delta = Delta_mat(i,j);
        numeric_eig = eig(jacobian_fcn(par.Delta,par.Gamma,par.Np,par.alpha,par.delta1,par.delta2,par.gammat));
        numeric_sol = sol_fcn(par.Delta,par.Gamma,par.Np,par.alpha,par.delta1,par.delta2,par.gammat);
        coex_stability(i,j) = (max(real((numeric_eig))) < 0) & (min(numeric_sol) > 0) & isreal(numeric_sol);
    end
end

predicted_trans_coex = total_trans_grid == 4;
agreement_transf = sum(sum(predicted_trans_coex == 0 & coex_stability == 0))/sum(sum(coex_stability==0));

%Create checkboard grid of states 1 and 2 to impose over the bistability
%(state 3)
bistable_pattern = double(checkerboard(50,10,10) > 0.5); 
bistable_pattern(bistable_pattern == 0) = 2;
total_trans_grid(total_trans_grid == 3) = bistable_pattern(total_trans_grid == 3);

%Plot
hold on
im2 = imagesc(total_trans_grid,[1,4]);
im2.AlphaData = par.alphaval;
plot(norm2im(Delta,Delta,par.n),norm2im(trans_plasmid_free_gamma,transfer_rates,par.n),'k-')
plot(norm2im(Delta,Delta,par.n),norm2im(trans_plasmid_only_gamma,transfer_rates,par.n),'k-')

xlabel({'$\log$(Plasmid cost, $\Delta$)'},'Interpreter','latex')
ylabel({'$\log$(Trans.','rate, $\gamma_\textrm{t}$)'},'Interpreter','latex')

xlim([1,par.n]);

ylim([1,par.n]);

yticks([])
xticks([])
box on
set(gca,'FontSize',4)

text(-0.25,1.12,labels{2},'Interpreter','latex','Units','normalized','FontSize',4)
end

