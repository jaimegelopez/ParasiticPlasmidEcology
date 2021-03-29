function pmf = plot_pmf(obs_vec,max_plasmids,max_n)

%This script plots a pmf

[counts,elements] = count_elements(obs_vec,max_plasmids);

pmf = counts/sum(counts);
markersize = 0.7;
hold on
plot(elements(1:max_n),pmf(1:max_n),'k-')
plot(elements(1:max_n),pmf(1:max_n),'ko','MarkerSize',markersize)

hold on
set(gca,'YScale','log')
xlim([-1,max_n])
xticks([0,max_n-1])
yticks([1e-4,1])
ylim([1e-4,1])
xticklabels({'0',num2str(max_n)})
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4);


end

