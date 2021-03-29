clear;clc

%This script generates figure 2

newfigure(3.42/2, 3.42/3*(2.3/2)*1.4);

gap = 0.1;
marg_h = [0.14,0.35];

Fig2ax = tight_subplot(2,2,[gap,0.2],marg_h,[0.15,0.05]);

fontsize = 4;

%% Leave space for biorender figure
axes(Fig2ax(1));

box off
text(-0.45,1.2,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)
axis off
axes(Fig2ax(2));
box off 
axis off


%% Plot different forms of epistasis

max_p = 9;
x = 0:(max_p-1);
Delta = 0.05;
axes(Fig2ax(3));

no_ep = 1- (1-Delta).^x; 
neg_ep = 1-(1-Delta).^(x.^1.5);
pos_ep = 1-(1-Delta).^(x > 0);

hold on
plot(x,no_ep,'g-')
plot(x,neg_ep,'b-')
plot(x,pos_ep,'r-')

ylim([0,0.8])
xlim([0,max_p])
xticks([0,4,8])
xticklabels({'0','4','8'})
yticks([0,0.4,0.8])
yticklabels({'0','0.4','0.8'})
xlabel({'Unique plasmid','number, $m$'},'Interpreter','latex','FontSize',fontsize)
ylabel('Total cost, $\Delta_{\textrm{tot}}$','Interpreter','latex','FontSize',fontsize)
text(-0.45,1.17,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')

leg_top = 0.77; 
leg_left = 0.4;
mylines = {'r-','g-','b-'};
labels = {'Pos. epis.','No epis.','Neg. epis.'};
line_length = 0.5;
spacing_y = 0.15;
spacing_x = 0.3;
curr_y = leg_top;
for i = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{i})
    text(leg_left + line_length + spacing_x,curr_y,labels{i},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.11;
end


%% Compute the different pmfs

plasmid_cost = 0.01;
n = 1000;
q = 0.005;
kappa = q/(1-q);

ii = 0:(n-1);

neutral_beta = (1-plasmid_cost).^ii;
neg_beta = (1-plasmid_cost).^(ii.^(1.5));

pos_beta = ones(1,n).*(1-plasmid_cost);
pos_beta(1) = 1;

neutral_pmf = compute_pmf(neutral_beta,kappa);
neg_pmf = compute_pmf(neg_beta,kappa);
pos_pmf = compute_pmf(pos_beta,kappa);

plot_max = max_p;
plot_n = ii(1:plot_max);

neutral_plot = neutral_pmf(1:plot_max);
neg_plot = neg_pmf(1:plot_max);
pos_plot = pos_pmf(1:plot_max);

markersize = 0.7;
axes(Fig2ax(4));
hold on
plot(plot_n,neutral_plot,'g-')
plot(plot_n,neutral_plot,'go','MarkerSize',markersize)

plot(plot_n,neg_plot,'b-')
plot(plot_n,neg_plot,'bo','MarkerSize',markersize)

plot(plot_n,pos_plot,'r-')
plot(plot_n,pos_plot,'ro','MarkerSize',markersize)

set(gca,'YScale','log')
ylim([0.5*min(neg_plot),1])
xlim([0,plot_max])
xticks([0,4,8])
xticklabels({'0','4','8'})
yticks([1e-8,1e-4,1e0])
yticklabels({'$10^{-8}$','$10^{-4}$','$10^0$'})
xlabel({'Unique plasmid','number, $m$'},'Interpreter','latex','FontSize',fontsize)
ylabel('Probability','Interpreter','latex','FontSize',fontsize)
text(-0.5,1.17,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca,'YMinorTick','Off')


print(gcf, '-dpng','fig2.png','-r1200');

close all