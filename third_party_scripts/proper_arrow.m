function [] = proper_arrow(x,y)

[figx1, figy1] = axxy2figxy(gca, x(1), y(1));
[figx2, figy2] = axxy2figxy(gca, x(2), y(2));


annotation('arrow',[figx1, figx2],[figy1, figy2],...
    'HeadLength',1.5,'HeadWidth',1.5)

end

