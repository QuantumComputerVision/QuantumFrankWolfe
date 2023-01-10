function plot_multiple(xValues, yMatrix, clrs, markers, dispNames, plotLegend)

for experiment=1:size(yMatrix,1)
    hold on, plot(xValues, yMatrix(experiment,:), 'LineWidth', 2, 'Color', clrs(experiment,:), 'DisplayName', dispNames{experiment}, 'Marker', markers{experiment});
end

if (plotLegend)
    lh=legend('show');
    set(lh,'FontSize',10);
end

end