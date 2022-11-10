% Plot a bar chart showing the proportions of photons ending in each final
% state

function [] = plot_probability_finalStates(final_state,inputs)


figure;
X = categorical({'Scatter out top','Scatter out Bottom','Absorbed'});
b = bar(X,[final_state.scatter_out_top, final_state.scatter_out_bottom, final_state.absorbed]);
set(gca,'TickLabelInterpreter', 'latex','FontSize',25);
set(gca, 'TitleFontSizeMultiplier',1.3)
bar_label_xLocation = b.XEndPoints;
bar_label_yLocation = b.YEndPoints;
bar_labels = string(b.YData./inputs.N_photons);
text(bar_label_xLocation,bar_label_yLocation,bar_labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',20,'FontWeight','bold','Interpreter','latex')
grid on; grid minor
ylabel('Counts','Interpreter','latex')
title('Probability of each final state', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1000 630])
ylim([0 inputs.N_photons])


dim = [0.7 0.85 0 0];
texBox_str = {['$N_{photons}^{total} = $', num2str(inputs.N_photons)],...
    ['$\lambda$ = ',num2str(inputs.wavelength), ' $nm$'],...
    ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
    ['$g$ = ', num2str(inputs.g)],...
    ['$r$ = ', num2str(inputs.radius), ' $\mu m$'],...
    ['$\tau_0$ = ', num2str(inputs.tau_upper_limit)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'white';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'white';
t.FitBoxToText = 'on';


end