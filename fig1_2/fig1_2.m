clear;clc

%This script generates the second part of figure 1

%% Simulate for figure

%Set-up parameters
par.n_plasmid = 8;
par.delta1 = 1;
par.f = 1e-1;
par.rl = 5e-2;
par.alpha1 = 1;
par.Gamma = 1;
par.error_threshold = 1e-6;
gammac_crit = (par.delta1/par.Gamma)*(par.delta1*par.f + par.delta1*par.rl*(1-par.f));
par.additive_fitness = 0;

%Variable conjugation rate
mult_vec =  [0.5, 2, 8, 32, 128];
gammac_vec = mult_vec.*gammac_crit;
count_cell1 = cell(size(gammac_vec));

for i = 1:length(gammac_vec)
    par.gammac = gammac_vec(i);
    par.dt = 0.1*min([1/par.alpha1,1/par.gammac]);
    count_cell1{i} = simulate_condensed_plasmid(par);
end

%Variable plasmid number
num_vec = [1 2 4 6 8];
max_num = max(num_vec);
par.gammac = 2*gammac_crit;    
par.dt = 0.1*min([1/par.alpha1,1/par.gammac]);
count_cell2 = cell(size(num_vec));

for i = 1:length(num_vec)
    par.n_plasmid = num_vec(i);
    count_cell2{i} = simulate_condensed_plasmid(par);
    count_cell2{i} = [count_cell2{i}; zeros(max_num-num_vec(i),1)];
end

save('fig1_2_data.mat')

%% Figure set-up

newfigure(3.42/2, 3.42/3*(2.4/2));

gap = [0.2,0.15];
marg_h = 0.1;

Fig1_2ax = tight_subplot(2,1,gap,marg_h,[0.14,0.05]);

fontsize = 4;

load('fig1_2_data.mat')

%% Make schematic

fig1_2_schematic(Fig1_2ax,1,{'\textbf{A}','\textbf{B}'},fontsize);


%% Variable conjugation rate histograms 

axes(Fig1_2ax(2));
hold on

%bar_color = [1,0,0];
bar_color = [0.9,0.9,0.9];
face_alpha = 0.3;

%Plot
space_mult = 1.4;
bar_labels1 =cellfun(@num2str,num2cell(0:par.n_plasmid),'UniformOutput',false);
tick_labels1 = cellfun(@num2str,num2cell(mult_vec),'UniformOutput',false);
tick_labels1 = strcat(tick_labels1,'$\gamma_\textrm{c}^*$');
base_x = 1:(par.n_plasmid + 1);
x_spacing = 0.05;
center_vec1 = zeros(size(gammac_vec));
for i = 1:length(gammac_vec)
    center_vec1(i) = median(base_x);
    bar(base_x,count_cell1{i},1,'FaceColor',bar_color,'FaceAlpha',face_alpha)
    for j = 1:length(bar_labels1)
       text(base_x(j),-0.05,bar_labels1{j},... %count_cell1{i}(j)+x_spacing
           'HorizontalAlignment','center','FontSize',2,'Interpreter','latex');
    end
    base_x = base_x + space_mult*par.n_plasmid;
end

yticks([0,1])
yticklabels({'0','1'})
xticks(center_vec1)
xticklabels(tick_labels1);
set(gca,'TickLength',[0,0]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out')
ylabel('Pop. fraction','Interpreter','latex')

text(-0.13,1.15,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)
text(-0.1,-0.17,'$\gamma_\textrm{c} = $','Interpreter','latex','Units','normalized','FontSize',4)


%% Save and close

print(gcf, '-dpng','fig1_2.png','-r1200','-painters');

close all



