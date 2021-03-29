function fig1_1_schematic(Fig1ax,loc,labels)

%% Plot the plasmid division event

axes(Fig1ax(loc(1)));

pos = get(gca,'Position');
pos(1) = pos(1) + pos(3)/1.5;
set(gca,'Position',pos);

plasmid_size = 0.05;
make_plasmid = @(coord)rectangle('Position',...
    [coord(1)-0.7*plasmid_size/2,coord(2)-1.3*plasmid_size/2,0.7*plasmid_size,1.3*plasmid_size],...
    'Curvature',[1,1],'LineWidth',0.5,'EdgeColor','r');

big_microbe_width = 2*0.25;
big_microbe_length = 2*0.07;
make_big_microbe = @(coord)rectangle('Position',...
    [coord(1)-big_microbe_length/2,coord(2)-big_microbe_width/2,big_microbe_length,big_microbe_width],...
    'Curvature',[1,0.5],'LineWidth',0.5);

coord_arrow = @(pt1,pt2) proper_arrow([pt1(1),pt2(1)],[pt1(2),pt2(2)]);
offset1 = [0.1,0];

big_offset = [-0.25 -0.15];
big1 = [0.25,0.5] + big_offset;
big1_1 = [0.75,0.8]+ big_offset;
big1_2 = [0.75,0.20]+ big_offset;

big2_1 = [1.25,0.8] + big_offset;
big2_2 = [1.25,0.2] + big_offset;

make_big_microbe(big1);
make_big_microbe(big1_1);
make_big_microbe(big1_2);
make_big_microbe(big2_1);
make_big_microbe(big2_2);

make_plasmid(big1);
make_plasmid(big1  + [0.04,0.06]);
make_plasmid(big1  + [-0.03,0.08]);
make_plasmid(big1  + [-0.01,0.17]);
make_plasmid(big1  + [-0.03,-0.17]);
make_plasmid(big1  + [0.04,-0.055]);
make_plasmid(big1  + [0.03,-0.15]);

make_plasmid(big1_1 + [0,0.1]);
make_plasmid(big1_1 + [-0.04,-0.1]);
make_plasmid(big1_1 + [0,-0.17]);

make_plasmid(big1_2  + [0.04,0.0]);
make_plasmid(big1_2  + [-0.03,0.04]);
make_plasmid(big1_2  + [-0.01,0.14]);
make_plasmid(big1_2  + [0.03,-0.17]);

make_plasmid(big2_1 + [-0.02,0.1]);
make_plasmid(big2_1 + [-0.04,-0.1]);
make_plasmid(big2_1 + [0,-0.17]);
make_plasmid(big2_1 + [0.04,-0.09]);
make_plasmid(big2_1 + [0.04,0.02]);
make_plasmid(big2_1 + [-0.02,0.01]);
make_plasmid(big2_1  + [0.03,0.17]);

make_plasmid(big2_2  + [0.04,0.0]);
make_plasmid(big2_2  + [-0.03,0.04]);
make_plasmid(big2_2  + [-0.03,0.14]);
make_plasmid(big2_2  + [0.03,-0.17]);
make_plasmid(big2_2  + [0.03,0.17]);
make_plasmid(big2_2  + [0,-0.04]);
make_plasmid(big2_2  + [-0.03,-0.15]);

xlim([0,1])
ylim([0,1])
axis off

coord_arrow(big1+offset1,big1_1 - offset1)
coord_arrow(big1+offset1,big1_2 - offset1)
coord_arrow(big1_1+offset1,big2_1 - offset1)
coord_arrow(big1_2+offset1,big2_2 - offset1)

ax = gca;    
ax.Clipping = 'off';   

text(-0.1 + big_offset(1),1.1,labels{1},'Interpreter','latex','Units','normalized','FontSize',4)

axes(Fig1ax(loc(2)));
axis off

%% Plot the conjugation schematic

axes(Fig1ax(loc(3)));

microbe_width = 1.4*0.25;
microbe_length = 1.4*0.07;
make_microbe = @(coord)rectangle('Position',...
    [coord(1)-microbe_length/2,coord(2)-microbe_width/2,microbe_length,microbe_width],...
    'Curvature',[1,0.5],'LineWidth',0.5);

xlim([0,1])
ylim([0,1])
axis off

schematic_spacing = 0.5;

p_offset = [0,plasmid_size];

primary = [0.2,0.3];
make_microbe(primary)
make_plasmid(primary+p_offset)
make_plasmid(primary-p_offset)

daughter_x_offset = 0.06;
microbe_spacing = 0.55;

horz_offset = [0,1.25*microbe_width/2];

death1 = [primary(1)+0.5*microbe_spacing,primary(2)];
text(death1(1),death1(2),'\boldmath$\emptyset$','HorizontalAlignment','center',...
    'VerticalAlignment','middle','Interpreter','latex','FontSize',4);

vert_offset = [1.4*microbe_length/2,0];
coord_arrow(primary+vert_offset,death1-vert_offset);

text(-0.1,1.1,labels{2},'Interpreter','latex','Units','normalized','FontSize',4)
text(0.5,1,'Conjugation','Interpreter','latex','Units','normalized',...
    'FontSize',4,'HorizontalAlignment','center','VerticalAlignment','middle')

empty_primary = [primary(1) + schematic_spacing,primary(2) + 0.25];
make_microbe(empty_primary);
donor_primary = [empty_primary(1)+2*daughter_x_offset,primary(2) + 0.25];
make_microbe(donor_primary);
make_plasmid(empty_primary+p_offset);
make_plasmid(empty_primary-p_offset);

transconjugant1 = [empty_primary(1),empty_primary(2)-1.2*microbe_spacing];
transconjugant2 =  [empty_primary(1)+2*daughter_x_offset,empty_primary(2)-1.2*microbe_spacing];
make_microbe(transconjugant1);
make_plasmid(transconjugant1+p_offset);
make_plasmid(transconjugant1-p_offset);

make_microbe(transconjugant2);
make_plasmid(transconjugant2+p_offset);
make_plasmid(transconjugant2-p_offset);

start1 = (empty_primary + donor_primary)/2;
end1 = (transconjugant1 + transconjugant2)/2;
coord_arrow(start1 - horz_offset,end1 + horz_offset);

ax = gca;    
ax.Clipping = 'off';   

%% Transformation schematic

axes(Fig1ax(loc(4)));

xlim([0,1])
ylim([0,1])
axis off

make_microbe(primary)
make_plasmid(primary+p_offset)
make_plasmid(primary-p_offset)


make_plasmid(death1+p_offset);
make_plasmid(death1-p_offset);

coord_arrow(primary+vert_offset,death1-vert_offset);

make_plasmid(empty_primary);
make_microbe(donor_primary);

make_microbe(end1);
make_plasmid(end1+p_offset);
make_plasmid(end1-p_offset);

coord_arrow(start1 - horz_offset,end1 + horz_offset);

text(-0.1,1.1,labels{3},'Interpreter','latex','Units','normalized','FontSize',4)

text(0.5,1,'Transformation','Interpreter','latex','Units','normalized',...
    'FontSize',4,'HorizontalAlignment','center','VerticalAlignment','middle')

ax = gca;    
ax.Clipping = 'off';   

end

