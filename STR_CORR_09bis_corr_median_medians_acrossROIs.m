function [corr_acrossROIs_WM_angle] = STR_CORR_09bis_corr_median_medians_acrossROIs(median_R1_WM,median_angle_WM,outpath,confidence_bounds)

% calculates correlation coefficient matrix computing the median of the
% medians across subjects and creates scatterplots across ROIs

N_subjs = size(median_R1_WM,1);

corr_acrossROIs_WM_angle = zeros(1,2);

plotname_WM = 'Scatterplot_WM_angle_vs_R1_acrossROIs';

median_median_R1_WM = median(median_R1_WM,1,'omitnan');
median_median_angle_WM = median(median_angle_WM,1,'omitnan');

mad_median_R1_WM = mad(median_R1_WM,1,1);
mad_median_angle_WM = mad(median_angle_WM,1,1);

SE_R1_WM = mad_median_R1_WM.*1.4826./sqrt(N_subjs);
SE_angle_WM = mad_median_angle_WM.*1.4826./sqrt(N_subjs);

h = figure;

% median across subjects for each region

[r,p] = corrcoef(median_median_R1_WM, median_median_angle_WM, 'rows','complete');

corr_acrossROIs_WM_angle(1) = r(2);
corr_acrossROIs_WM_angle(2) = p(2);

corr_coef = round(r(2),3,'significant');
p_value = round(p(2),3,'significant');

corr_coef_str_WM = num2str(corr_coef);
p_value_str_WM = num2str(p_value);

h.Color = [1 1 1];
h.Position = [1 1 1200 800];

s = errorbar(median_median_R1_WM, median_median_angle_WM, SE_angle_WM, SE_angle_WM, SE_R1_WM, SE_R1_WM,'o');
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'r';
s.MarkerFaceColor = [0.8 0.8 0.8];
s.Color = 'r';
hold on
plot(median_median_R1_WM, median_median_angle_WM,'ro');


if confidence_bounds == 1
    % confidence bounds
    x = median_median_R1_WM;
    y = median_median_angle_WM;
    [p,S] = polyfit(x,y,1);
    % Ensure the confidence interval covers the full x-axis
    xLimits = xlim;  % Get current x-axis limits
    xv = linspace(xLimits(1), xLimits(2), 150);  % Generate x values over full x-axis
    [y_ext,delta] = polyconf(p,xv,S);
    hold on
    %plot(xv, y_ext, '-b')
    patch([xv fliplr(xv)], [(y_ext+delta) fliplr((y_ext-delta))], 'b', 'FaceAlpha',0.1, 'EdgeColor','none')
end

hold on

kk = lsline;
kk.Color = 'r'; kk.LineWidth = 2;
ylabel("Angle (°)",'FontSize',14,'FontWeight','Bold');
xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
%xlim([0.5 1]); ylim([0 0.5]);

title({strcat('r_{WM}=',corr_coef_str_WM,'; p=',p_value_str_WM)});

ax = gca;
ax.Box = 'off'; % Remove the box around the plot
ax.FontSize = 12;


saveas(gcf,fullfile(outpath,plotname_WM),'fig');


end