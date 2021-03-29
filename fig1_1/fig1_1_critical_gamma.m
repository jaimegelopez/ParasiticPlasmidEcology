function fig1_1_critical_gamma(Fig1ax,par,loc,labels)

%% Set-up

pcns = linspace(par.min_pcn,par.max_pcn,par.n);

Delta = pcns.*par.per_plasmid_cost;
rl = 2.^(1-pcns);

conj_fit_loss = par.delta1.*Delta;
conj_seg_loss = par.delta1*(1 - Delta).*rl;

trans_denom = pcns - Delta - rl.*(1-Delta);
trans_fit_loss = par.delta2.*Delta./trans_denom; 
trans_seg_loss = par.delta2.*rl.*(1-Delta)./trans_denom;

color = [1,0,0];
FaceAlpha = 0.3;
%% Plot conjugation loss
axes(Fig1ax(loc(1)));
hold on
conj_total_loss = (par.delta1/par.Gamma)*(conj_fit_loss + conj_seg_loss);
plot(pcns,conj_total_loss,'k-')

ylim([0,par.loss_ymax])
xlim([par.min_pcn,par.max_pcn])
xticks([])
yticks([])
xlabel('Copy number, $n_\textrm{p}$','Interpreter','latex')
ylabel({'Conjug.', 'rate, $\gamma_\textrm{c}$'},'Interpreter','latex')
set(gca,'FontSize',4)
text(-0.25,1.12,labels{1},'Interpreter','latex','Units','normalized','FontSize',4)
box on

patch([pcns fliplr(pcns)], [conj_total_loss par.loss_ymax*ones(size(pcns))],...
    color,'FaceAlpha',FaceAlpha)        

text(0.35,0.7,{'Plasmid','invasion'},'Interpreter','latex','Units','normalized',...
    'FontSize',4,'HorizontalAlignment','center','VerticalAlignment','middle')


%% Plot transformation loss
axes(Fig1ax(loc(2)));
hold on
trans_total_loss = (par.delta1/par.Gamma)*(trans_fit_loss + trans_seg_loss);
plot(pcns,trans_total_loss,'k-')

ylim([0,par.loss_ymax])
xlim([par.min_pcn,par.max_pcn])
xticks([])
yticks([])
xlabel('Copy number, $n_\textrm{p}$','Interpreter','latex')
ylabel({'Trans.', 'rate, $\gamma_\textrm{t}$'},'Interpreter','latex')
set(gca,'FontSize',4)

text(-0.25,1.12,labels{2},'Interpreter','latex','Units','normalized','FontSize',4)
patch([pcns fliplr(pcns)], [trans_total_loss par.loss_ymax*ones(size(pcns))],...
    color,'FaceAlpha',FaceAlpha)        

box on

text(0.5,0.45,{'Plasmid','invasion'},'Interpreter','latex','Units','normalized',...
    'FontSize',4,'HorizontalAlignment','center','VerticalAlignment','middle')


end

