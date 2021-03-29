clear;clc

%This script removes major genera and performs the distribution analysis

load('fig3_analysis.mat');

count_cutoff = 10; %Truncate pmfs when there are <10 observations

%Manually curated list of genera
clinical_genera = {'escherichia','klebsiella','shigella','chlamydia',...
    'streptococcus','pseudomonas','staphylococcus','micrococcus','bacillus',...
    'corynebacterium','listeria','nocardia',...
    'gardnerella','lactobacillus','actinomyces','clostridium','propionibacterium',...
    'moraxella','neisseria','veillonella','haemophilus','vibrio','salmonella',...
    'acinetobacter','campylobacter','burkholderia','aeromonas','proteus',...
    'stenotrophomonas','bacteroides','fusobacterium','prevotella','enterobacter',...
    'citrobacter','enterococcus','listeria','leptospira',...
    'mycoplasma','pasteurella','bordetella','borrelia','yersinia','legionella',...
    'mycobacterium','enterobacteriaceae','helicobacter','brucella','serratia',...
    'francisella','rickettsia','clostridioides','chryseobacterium','flavobacterium',...
    'mycobacteroides','elizabethkingia','microbacterium','pantoea','bartonella',...
    'sphingomonas','arcobacter','providencia','achromobacter','borreliella',...
    'capnocytophaga','pandoraea','histophilus','leclercia','actinobacillus',...
    'treponema','cronobacter','coxiella','edwardsiella','glaesserella',...
    'ureaplasma','candidatus','bifidobacterium','akkermansia','lactococcus',...
    'leuconostoc','cutibacterium','desulfovibrio','alistipes','lachnospiraceae',...
    'myroides','rothia','blautia','faecalibacterium','eubacterium','hungateiclostridium',...
    'eggerthella','collinsella','alcaligenes'}';
clinical_genera = unique(clinical_genera);

%Loop through and remove genera
total_ind = zeros(size(filt_origin_table,1),1);
for i = 1:length(clinical_genera)
    genera_i = clinical_genera{i};
    ind_i = strcmp(filt_origin_table.genus,genera_i)|strcmp(filt_origin_table.genus,['[',genera_i,']']);
    total_ind = total_ind | ind_i;
end
nonpath_table = filt_origin_table(~total_ind,:);
[genus,genus_counts] = count_taxa_frequency(nonpath_table,'genus');

%Compute the all genome pmf from the data
alphaval = 1;
fontsize = 10;
[nonpath_pmf,xvec] = get_data_pmf(nonpath_table.n_plasmids,count_cutoff);
nonpath_max = length(nonpath_pmf)-1;
x = (0:nonpath_max)';
%Set up exponent vectors for fitting
n = 1000;
pos_exp_vec = [0;ones(n-1,1)];
[pos_pmf, pos_par] = fit_pmf_to_data(nonpath_pmf,pos_exp_vec);

%Plot the supp figure
newfigure(3.42/2, 3.42/3*(2/2));
hold on
patchline(x,nonpath_pmf,'linestyle','-','edgecolor','k');
scatter(x,nonpath_pmf,0.7,'ko');
patchline(x,pos_pmf(x+1),'linestyle','-','edgecolor','r','edgealpha',alphaval);
scatter(x,pos_pmf(x+1),0.7,'o','MarkerEdgeAlpha',alphaval,'MarkerFaceColor','r','MarkerEdgeColor','r');
text(0.05,0.15,'Non-clinical spp.','Units','normalized','FontSize',4,'Interpreter','latex')
xlabel({'Unique plasmid number, $m$'},'Interpreter','latex','FontSize',fontsize)
ylabel({'Probability'},'Interpreter','latex','FontSize',fontsize)
set(gca,'YScale','log')
xlim([-1,length(x)+1])
xticks([0,5,length(x)+1])
genus_ylim = [1e-3,1e0];
genus_ytick = [1e-3,1e-2, 1e-1,1e0];
genus_ylabel = {'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^0$'};
yticks(genus_ytick)
ylim(genus_ylim)
xticklabels({'0','5','10'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4);
set(gca,'YMinorTick','Off')

print(gcf, '-dpng','fig3_nonclinical.png','-r1200');
