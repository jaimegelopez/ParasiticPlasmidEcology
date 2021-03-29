clear;clc
load('fig3_analysis.mat');

%% Set-up
count_cutoff = 10; %Truncate pmfs when there are <10 observations

%Set up exponent vectors for fitting
n = 1000;
neut_exp_vec = (0:n-1)';
alphaval = 1;

%set-up shared plot properties
genus_ytick = [1e-2,1e-1,1e0];
genus_ylabel = {'$10^{-2}$','$10^{-1}$','$10^0$'};

%Get top genera we didn't cover in fig 4
large_genera = genus([3:7,9:12]);

labels = {'A','B','C','D','E','F','G','H','I'};

%% Loop through the genera and fit models

for i = 1:length(large_genera)
    genus = large_genera{i};
    taxa_table = filt_origin_table(strcmp(filt_origin_table.genus,genus),:);
    [gen_pmf{i},gen_x{i}] = get_data_pmf(taxa_table.n_plasmids,count_cutoff);
    gen_max(i) = length(gen_pmf{i}) - 1;
    [pos_gen_pmf{i},pos_gen_par{i}] = fit_var_epis_pmf(gen_pmf{i},neut_exp_vec);
    
end

save('genus_level_pmfs.mat')

%% Loop through and plot

load('genus_level_pmfs.mat')

newfigure(3.42/2*1.5, 3.42/3*(3.5/2));

gap = [0.13,0.05];
marg_h = [0.15,0.11];

supp_Fig4ax = tight_subplot(3,3,gap,marg_h,[0.13,0.13]);

fontsize = 4;

convstr = @(x) num2str(round(x,2,'significant'));

for i = 1:length(supp_Fig4ax)
    axes(supp_Fig4ax(i));
    hold on
    %Plot genus pmfs
    genus = large_genera{i};
    patchline(gen_x{i},gen_pmf{i},'linestyle','-','edgecolor','k');
    scatter(gen_x{i},gen_pmf{i},0.7,'ko');
    patchline(gen_x{i},pos_gen_pmf{i}(gen_x{i}+1),'linestyle','-','edgecolor','r','edgealpha',alphaval);
    scatter(gen_x{i},pos_gen_pmf{i}(gen_x{i}+1),0.7,'o','MarkerEdgeAlpha',alphaval,'MarkerFaceColor','r','MarkerEdgeColor','r');
    gen_label = genus;
    gen_label(1) = upper(gen_label(1));
    text(0.05,1.1,['\textit{',gen_label,'} spp.'],'Units','normalized','FontSize',4,'Interpreter','latex')
    
    yticklabels(genus_ylabel)
    text(-0.15,1.25,['\textbf{',labels{i},'}'],'Interpreter','latex','Units','normalized','FontSize',4)
    set(gca,'YScale','log')
    xlim([-1,gen_max(i)+1])
    xlim([-1,9])
    ylim([1e-2,1e0])
    yticks(genus_ytick)
    xticks([0,4,8])
    
    if sum(i == [1,4,7])
        ylabel({'Probability'},'Interpreter','latex','FontSize',fontsize)
    else
        yticklabels([])
    end
    if sum(i == [7,8,9])
        xlabel({'Unique plasmid','number, $m$'},'Interpreter','latex','FontSize',fontsize)
        xticklabels({'0','4','8'})
    end
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',4);
    set(gca,'YMinorTick','Off')
    
    if i == 3
        %Generate legend
        leg_top = 5e0;
        leg_left = 9;
        mylines = {'k-','r-'};
        legend_labels = {'Data','Model fit'};
        line_length = 0.5;
        spacing_x = 0.3;
        curr_y = leg_top;
        for j = 1:(length(legend_labels))
            plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{j})
            text(leg_left + line_length + spacing_x,curr_y,legend_labels{j},'Interpreter','latex','FontSize',4)
            curr_y = curr_y*5e-1;
        end
        ax = gca;
        ax.Clipping = 'off';
    end
    
        leg_top = 0.4e0;
        leg_left = 4;
        legend_labels = {['$\Delta = ',convstr(pos_gen_par{i}(1)),'$'],...
            ['$q = ',convstr(pos_gen_par{i}(2)),'$'],...
            ['$a = ',convstr(pos_gen_par{i}(3)),'$']};
        line_length = 0.5;
        spacing_x = 0.3;
        curr_y = leg_top;
        for j = 1:(length(legend_labels))
            text(leg_left + line_length + spacing_x,curr_y,legend_labels{j},'Interpreter','latex','FontSize',4)
            curr_y = curr_y*5e-1;
        end
        ax = gca;
        ax.Clipping = 'off';
    
end

print(gcf, '-dpng','supp_figure_genus_level_pmfs.png','-r1200');
