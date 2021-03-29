function fig2_schematic(Fig2ax,loc,labels,fontsize)

%This function makes the multiplasmid schematic

axes(Fig2ax(loc(1)));

ora = 'b';

plasmid_size = 0.05;
make_plasmid = @(coord,color)rectangle('Position',...
    [coord(1)-0.7*plasmid_size/2,coord(2)-1.3*plasmid_size/2,0.7*plasmid_size,1.3*plasmid_size],...
    'Curvature',[1,1],'LineWidth',0.5,'EdgeColor',color);

big_microbe_width = 2*0.25;
big_microbe_length = 2*0.07;
make_big_microbe = @(coord)rectangle('Position',...
    [coord(1)-big_microbe_length/2,coord(2)-big_microbe_width/2,big_microbe_length,big_microbe_width],...
    'Curvature',[1,0.5],'LineWidth',0.5);

coord_arrow = @(pt1,pt2) proper_arrow([pt1(1),pt2(1)],[pt1(2),pt2(2)]);
offset1 = [0.1,0];

big_offset = [-0.4,0];
big1 = [0.25,0.5] + big_offset;
big1_1 = [0.65,0.8]+ big_offset;
big1_2 = [0.65,0.20]+ big_offset;
big2_1 = big1_1 + [0.4,0];
big2_2 = big1_2 + [0.4,0];

make_big_microbe(big1);
make_big_microbe(big1_1);
make_big_microbe(big1_2);
make_big_microbe(big2_1);
make_big_microbe(big2_2);

make_plasmid(big1,'r');
make_plasmid(big1  - [0.04,0.06],'r');
make_plasmid(big1  - [-0.03,0.08],ora);
make_plasmid(big1  - [-0.01,0.17],'r');
make_plasmid(big1  - [-0.03,-0.13],ora);
make_plasmid(big1  - [0.04,-0.055],'r');
make_plasmid(big1  - [0.03,-0.15],ora);

make_plasmid(big1_1 - [0,0.1],'r');
make_plasmid(big1_1 - [-0.04,-0.07],ora);
make_plasmid(big1_1 - [0,-0.14],'r');

make_plasmid(big1_2  - [0.04,0.0],ora);
make_plasmid(big1_2  - [-0.03,-0.04],'r');
make_plasmid(big1_2  - [-0.01,0.14],ora);
make_plasmid(big1_2  - [-0.03,-0.17],'r');

make_plasmid(big2_1 - [0,0.09],'r');
make_plasmid(big2_1 - [-0.035,0.14],ora);
make_plasmid(big2_1 - [0.035,0.14],'r');
make_plasmid(big2_1 - [-0.0,-0.0],ora);
make_plasmid(big2_1 - [-0.04,-0.07],ora);
make_plasmid(big2_1 - [0,-0.14],'r');
make_plasmid(big2_1 - [0.04,-0.07],'r');

make_plasmid(big2_2  - [0.04,0.0],ora);
make_plasmid(big2_2  - [-0.03,-0.04],'r');
make_plasmid(big2_2  - [-0.01,0.14],ora);
make_plasmid(big2_2  - [-0.03,-0.17],'r');
make_plasmid(big2_2  - [0.02,-0.16],ora);
make_plasmid(big2_2  - [0.02,-0.07],'r');
make_plasmid(big2_2  - [0.0,0.05],'r');


xlim([0,2])
ylim([0,1])
axis off

coord_arrow(big1+offset1,big1_1 - offset1)
coord_arrow(big1+offset1,big1_2 - offset1)
coord_arrow(big1_1+offset1,big2_1 - offset1)
coord_arrow(big1_2+offset1,big2_2 - offset1)

text(-0.13,1.12,labels{1},'Interpreter','latex','Units','normalized','FontSize',fontsize)

ax = gca;    
ax.Clipping = 'off';   

%% Make the plasmid cost figure


pos = get(gca,'Position');
pos(1) = pos(2) + 0.01*pos(3);
pos(3) = 0.4*pos(3);

subax = axes('Position',pos);
xticks([])
yticks([])
set(gca,'FontSize',fontsize)

delta = 0.05;
xrange = 0:8;
ydata = (1 - delta).^xrange;

plot(xrange,1-ydata,'k-')

xlim([0,9])
ylim([0,1.1*max(1-ydata)])
xticks([])
yticks([])
box off
xlabel('Unique plasmid number, $m$','Interpreter','latex','FontSize',fontsize)
ylabel('Total cost, $\Delta_\textrm{tot}$','Interpreter','latex','FontSize',fontsize)
text(-0.35,1.12,labels{2},'Interpreter','latex','Units','normalized','FontSize',fontsize)

xticks([0,4,8])
xticklabels({'0','4','8'})
yticks([0,0.2,0.4])
yticklabels({'0','0.2','0.4'})
ylim([0,0.4])

set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')

ax = gca;    
ax.Clipping = 'off';   

end
