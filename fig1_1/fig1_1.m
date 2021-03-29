 clear;clc

%This script makes the first part of figure 1. This previously was a
%standalone figure 1, but is now merged into a combined figure 1

newfigure(3.42/2, (4/3)*3.42/2);

gap = [0.1,0.15];
marg_h = [0.13,0.04];
marg_w = [0.12,0.1];

Fig1_1ax = tight_subplot(4,2,gap,marg_h,marg_w);

par.n = 1000;

%% Plot the conjugation schematic

schem_labels = {'\textbf{A}','\textbf{B}','\textbf{E}'};
fig1_1_schematic(Fig1_1ax,[1,2,3,4],schem_labels);

%% Plot plasmid loss stuff

par.max_pcn = 15;
par.min_pcn = 1.1;

par.Gamma = 1;
par.delta1 = 1;
par.delta2 = 1;
par.alpha = 1;
par.loss_ymax = 1.1;

par.per_plasmid_cost = 0.05;

crit_gamma_labels = {'\textbf{C}','\textbf{F}'};

fig1_1_critical_gamma(Fig1_1ax,par,[5,6],crit_gamma_labels);


%% Plot phase diagrams

%Parameters
par.delta1 = 0.1;
par.alpha = 1;
par.delta2 = 0.3;
par.Np = 0.6; 

par.alphaval = 0.8;

phase_labels = {'\textbf{D}','\textbf{G}'};

[agreement_conj,agreement_transf] = fig1_1_phase_diagram(Fig1_1ax,par,[7,8],phase_labels);


%% Print file

print(gcf, '-dpng','fig1_1.png','-r1200','-painters');

close all