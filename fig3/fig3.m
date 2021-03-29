clear;clc

%This script generates figure 3

%% Read info

filter_bioproject = false;

element_info = readtable('merged_element_info.csv','Delimiter','\t');

element_info = fix_gene_lists(element_info,{'cas'});

%% Get the origin table

origin_table = get_origin_table(element_info);

%% Filter based on Bioproject id

if filter_bioproject
    filt_origin_table = filter_on_bioproject(origin_table);
else
    filt_origin_table = origin_table;
end

%% Get genera and species

filt_origin_table = assign_taxonomy(filt_origin_table);

%% Count species frequency

[species,species_counts] = count_taxa_frequency(filt_origin_table,'species');
taxa_cutoff = 250;
common_species = species(species_counts > taxa_cutoff);

[genus,genus_counts] = count_taxa_frequency(filt_origin_table,'genus');
common_genus = genus(genus_counts > taxa_cutoff);

save('fig3_analysis.mat');


%% Prepare plot
clear;clc
load('fig3_analysis.mat');
filt_origin_table.cas_present = logical(filt_origin_table.cas_present);

newfigure(3.42/2, 3.42/3*(3.5/2));

gap = [0.2,0.05];
marg_h = [0.1,0.07];

Fig4ax = tight_subplot(3,2,gap,marg_h,[0.16,0.05]);

fontsize = 4;

axes(Fig4ax(1));
pos = get(gca,'Position');
pos(1) = 1.3*pos(1);
pos(3) = 1.8*pos(3);
set(gca,'Position',pos);

axes(Fig4ax(5));
pos = get(gca,'Position');
pos(1) = 1.3*pos(1);
pos(3) = 1.8*pos(3);
pos(2)= 1.08*pos(2);
set(gca,'Position',pos);

axes(Fig4ax(2));
axis off

axes(Fig4ax(6));
axis off

%% Plot all plasmid distribution

axes(Fig4ax(1));
hold on
count_cutoff = 10; %Truncate pmfs when there are <10 observations

%Compute the all genome pmf from the data
[all_pmf,xvec] = get_data_pmf(filt_origin_table.n_plasmids,count_cutoff);
all_max = length(all_pmf)-1;

%Set up exponent vectors for fitting
n = 1000;
pos_exp_vec = [0;ones(n-1,1)];
neut_exp_vec = (0:n-1)';
neg_exp_vec = (0:n-1)'.^1.5;

%Fit positive and neutral epistasis fits to all genomes data
[pos_pmf, pos_par] = fit_pmf_to_data(all_pmf,pos_exp_vec);
[neut_pmf, neut_par] = fit_pmf_to_data(all_pmf,neut_exp_vec);
indvec = xvec + 1;
alphaval = 1;

%Plot fig A - all genomes and fits
p(1) = patchline(xvec,all_pmf,'linestyle','-','edgecolor','k');
scatter(xvec,all_pmf,0.7,'ko');
p(2) = patchline(xvec,pos_pmf(xvec+1),'linestyle','-','edgecolor','r','edgealpha',alphaval);
scatter(xvec,pos_pmf(xvec+1),0.7,'ro','MarkerEdgeAlpha',alphaval,'MarkerFaceColor','r','MarkerFaceAlpha',alphaval);
p(3) = patchline(xvec,neut_pmf(xvec+1),'linestyle','-','edgecolor','g','edgealpha',alphaval);
scatter(xvec,neut_pmf(xvec+1),0.7,'go','MarkerEdgeAlpha',alphaval,'MarkerFaceColor','g','MarkerFaceAlpha',alphaval);
set(gca,'YScale','log')
xlim([-1,all_max+1])
xticks([0,4,8,all_max+1])
yticks([1e-4,1e-2,1])
yticklabels({'$10^{-4}$','$10^{-2}$','$10^0$'})
ylim([1e-4,1])
xticklabels({'0','4','8',num2str(all_max+1)})
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4);
box off
xlabel({'Unique plasmid number, $m$'},'Interpreter','latex','FontSize',fontsize)
text(0.1,0.2,'All genomes','Units','normalized','FontSize',4,'Interpreter','latex')
ylabel({'Probability'},'Interpreter','latex','FontSize',fontsize)
text(-0.2,1.25,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'YMinorTick','Off')

%Generate legend
leg_top = 5e0; 
leg_left = 7.5;
mylines = {'k-','r-','g-','b-'};
labels = {'Data','Pos. epistasis','No epistasis','Neg. epistasis'};
line_length = 0.5;
spacing_x = 0.3;
curr_y = leg_top;
for i = 1:(length(labels)-1)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{i})
    text(leg_left + line_length + spacing_x,curr_y,labels{i},'Interpreter','latex','FontSize',4)
    curr_y = curr_y*1.5e-1;
end
ax = gca;    
ax.Clipping = 'off';   


%% fig D - cas and non cas
gray = [0.5 0.5 0.5];

%Generate pmfs of cas containing and non-cas containing genomes
[cas_pmf,cas_x] = get_data_pmf(filt_origin_table.n_plasmids(filt_origin_table.cas_present),count_cutoff);
[non_cas_pmf, non_cas_x] = get_data_pmf(filt_origin_table.n_plasmids(~filt_origin_table.cas_present),count_cutoff);

%plot cas vs. non-cas
axes(Fig4ax(5));
hold on
patchline(cas_x,cas_pmf,'linestyle','-','edgecolor',gray);
scatter(cas_x,cas_pmf,0.7,'o','MarkerFaceColor',gray,'MarkerEdgeColor',gray);
patchline(non_cas_x,non_cas_pmf,'linestyle','--','edgecolor','k');
scatter(non_cas_x,non_cas_pmf,0.7,'ko');
set(gca,'YScale','log')
xlim([-1,length(cas_x)])
xticks([0,4,8,length(cas_x)])
yticks([1e-4,1e-2,1])
ylim([1e-4,1])
xticklabels({'0','4','8',length(cas_x)})
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4);
box off
xlabel({'Unique plasmid number, $m$'},'Interpreter','latex','FontSize',fontsize)
ylabel({'Probability'},'Interpreter','latex','FontSize',fontsize)
yticklabels({'$10^{-4}$','$10^{-2}$','$10^0$'})
text(-0.2,1.25,'\textbf{D}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'YMinorTick','Off')

%cas vs. non-cas legend
leg_top = 2e-1; 
leg_left = 8;
mylines = {gray,'k'};
linetypes = {'-','-'};
labels = {'\textit{cas}','No \textit{cas}'};
line_length = 0.5;
spacing_y = 0.15;
spacing_x = 0.3;
curr_y = leg_top;
for i = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],linetypes{i},'Color',mylines{i})
    text(leg_left + line_length + spacing_x,curr_y,labels{i},'Interpreter','latex','FontSize',4)
    curr_y = curr_y*2e-1;
end


%% fig B-C genus level 

%set-up shared plot properties
genus_ylim = [1e-2,1e0];
genus_ytick = [1e-2,1e-1,1e0];
genus_ylabel = {'$10^{-2}$','$10^{-1}$','$10^0$'};
axes(Fig4ax(3));
hold on

%Get escherichia genomes and compute pmf with arbitrary exponent
taxa_table = filt_origin_table(strcmp(filt_origin_table.genus,'escherichia'),:);
[ecoli_pmf,ecoli_x] = get_data_pmf(taxa_table.n_plasmids,count_cutoff);
ecoli_max = length(ecoli_pmf) - 1;
[pos_ecoli_pmf,pos_ecoli_par] = fit_var_epis_pmf(ecoli_pmf,neut_exp_vec);

%% Plot escherichia pmfs
patchline(ecoli_x,ecoli_pmf,'linestyle','-','edgecolor','k');
scatter(ecoli_x,ecoli_pmf,0.7,'ko');
patchline(ecoli_x,pos_ecoli_pmf(ecoli_x+1),'linestyle','-','edgecolor','r','edgealpha',alphaval);
scatter(ecoli_x,pos_ecoli_pmf(ecoli_x+1),0.7,'o','MarkerEdgeAlpha',alphaval,'MarkerFaceColor','r','MarkerEdgeColor','r');
text(0.05,0.15,'\textit{Escherichia} spp.','Units','normalized','FontSize',4,'Interpreter','latex')
xlabel({'Unique plasmid','number, $m$'},'Interpreter','latex','FontSize',fontsize)
ylabel({'Probability'},'Interpreter','latex','FontSize',fontsize)
yticklabels(genus_ylabel)
text(-0.35,1.25,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'YScale','log')
xlim([-1,ecoli_max+1])
xticks([0,4,ecoli_max+1])
yticks(genus_ytick)
ylim(genus_ylim)
xticklabels({'0','4',num2str(ecoli_max+1)})
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4);
set(gca,'YMinorTick','Off')

%% Get kleb genomes and fit pmf
taxa_table = filt_origin_table(strcmp(filt_origin_table.genus,'klebsiella'),:);
[kleb_pmf,kleb_x] = get_data_pmf(taxa_table.n_plasmids,count_cutoff);
kleb_max = length(kleb_pmf)-1;

%Plot kleb distribution
axes(Fig4ax(4));
hold on
[pos_kleb_pmf,pos_kleb_par] = fit_var_epis_pmf(kleb_pmf,neut_exp_vec);
patchline(kleb_x,kleb_pmf,'linestyle','-','edgecolor','k');
scatter(kleb_x,kleb_pmf,0.7,'ko');
patchline(kleb_x,pos_kleb_pmf(kleb_x+1),'linestyle','-','edgecolor','r');
scatter(kleb_x,pos_kleb_pmf(kleb_x+1),0.7,'o','MarkerFaceColor','r','MarkerEdgeColor','r');
text(0.05,0.15,'\textit{Klebsiella} spp.','Units','normalized','FontSize',4,'Interpreter','latex')
xlabel({'Unique plasmid','number, $m$'},'Interpreter','latex','FontSize',fontsize)
text(-0.15,1.25,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'YScale','log')
xlim([-1,kleb_max+1])
xticks([0,4,kleb_max+1])
yticks(genus_ytick)
ylim(genus_ylim)
xticklabels({'0','4',num2str(kleb_max+1)})
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4);
set(gca,'YMinorTick','Off')

if filter_bioproject
    print(gcf, '-dpng','fig3_supp/fig3_no_artificial.png','-r1200');
else
    print(gcf, '-dpng','fig3.png','-r1200');
end


