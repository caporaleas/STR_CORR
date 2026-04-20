function [corr_acrossROIs_noout_GM,corr_acrossROIs_noout_WM] = STR_CORR_11bis_corr_median_medians_acrossROIs_noout_IQR(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,outpath,confidence_bounds)

% calculates correlation coefficient matrix computing the median of the
% medians across subjects and creates scatterplots, after removing outliers
% for each parameter
N_par = length(parameters);
N_GM_ROIs = size(median_R1_GM,2);
N_WM_ROIs = size(median_R1_WM,2);
N_subjs = size(median_R1_GM,1);

corr_acrossROIs_noout_GM = zeros(N_par-1,2);
corr_acrossROIs_noout_WM = zeros(N_par-1,2);

plotname_GM = 'Scatterplot_GM_diffpar_vs_R1_acrossROIs_noout_IQR';
plotname_WM = 'Scatterplot_WM_diffpar_vs_R1_acrossROIs_noout_IQR';

median_median_R1_GM = median(median_R1_GM,1,'omitnan');
median_median_diff_GM = median(median_diff_GM,1,'omitnan');

mad_median_R1_GM = mad(median_R1_GM,1,1);
mad_median_diff_GM = mad(median_diff_GM,1,1);

median_median_R1_WM = median(median_R1_WM,1,'omitnan');
median_median_diff_WM = median(median_diff_WM,1,'omitnan');

mad_median_R1_WM = mad(median_R1_WM,1,1);
mad_median_diff_WM = mad(median_diff_WM,1,1);

for pp = 1 : N_par-1
    tmp = median_median_diff_GM(1,:, pp);
    out_idx = isoutlier(tmp,'quartiles');
    median_median_diff_GM(1,out_idx,pp) = NaN;
end

for pp = 1 : N_par-1
    tmp = median_median_diff_WM(1,:, pp);
    out_idx = isoutlier(tmp,'quartiles');
    median_median_diff_WM(1,out_idx,pp) = NaN;
end

SE_R1_GM = mad_median_R1_GM.*1.4826./sqrt(N_subjs);
SE_diff_GM = mad_median_diff_GM.*1.4826./sqrt(N_subjs);

SE_R1_WM = mad_median_R1_WM.*1.4826./sqrt(N_subjs);
SE_diff_WM = mad_median_diff_WM.*1.4826./sqrt(N_subjs);

correlations = zeros(N_par-1,1);
rawpvalues = zeros(N_par-1,1);

h = figure;
for pp = 1 : N_par-1
    
    subplot(2,4,pp) 

    % median across subjects for each region
       
    [r,p] = corrcoef(median_median_R1_GM, median_median_diff_GM(1,:,pp), 'rows','complete');

    correlations(pp,1) = r(2);
    rawpvalues(pp,1) = p(2);

    corr_acrossROIs_noout_GM(pp,1) = r(2);
    corr_acrossROIs_noout_GM(pp,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];
    h.Position = [1 1 1200 800];

    s = errorbar(median_median_R1_GM, median_median_diff_GM(1,:,pp), SE_diff_GM(1,:,pp), SE_diff_GM(1,:,pp), SE_R1_GM, SE_R1_GM,'o');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = [0.8 0.8 0.8];
    s.Color = 'b';
    hold on
    plot(median_median_R1_GM, median_median_diff_GM(1,:,pp),'bo');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_GM;
        y = median_median_diff_GM(1,:,pp);
        idx = isnan(y);
        y(idx) = [];
        x(idx) = [];
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
    kk.Color = 'b'; kk.LineWidth = 2;
    ylabel(parameters{1,pp+1},'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;

end
saveas(gcf,fullfile(outpath,plotname_GM),'fig');

% FDR correction using Benjamini-Hochberg method
if exist('mafdr', 'file') % Check if mafdr is available
    adjpvalues = mafdr(rawpvalues, 'BHFDR', true);
else
    % Manual FDR correction if mafdr is not available
    [sortedP, sortIdx] = sort(rawpvalues);
    numTests = length(pValues);
    adjustedP = sortedP .* numTests ./ (1:numTests)';
    adjustedP(adjustedP > 1) = 1; % Ensure p-values are not >1
    
    % Reorder back to original indices
    adjpvalues = zeros(size(rawpvalues));
    adjpvalues(sortIdx) = adjustedP;
end

length(correlations)
resultsTable = table(parameters(2:end)',correlations, rawpvalues, adjpvalues, ...
                     'VariableNames', {'Parameter','GM_Correlation', 'Raw_pValue', 'Adjusted_pValue'});
disp(resultsTable);


correlations = zeros(N_par-1,1);
rawpvalues = zeros(N_par-1,1);

h = figure;
for pp = 1 : N_par-1
    
    subplot(2,4,pp) 

    % median across subjects for each region
       
    [r,p] = corrcoef(median_median_R1_WM, median_median_diff_WM(1,:,pp), 'rows','complete');

    correlations(pp,1) = r(2);
    rawpvalues(pp,1) = p(2);

    corr_acrossROIs_noout_WM(pp,1) = r(2);
    corr_acrossROIs_noout_WM(pp,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_WM = num2str(corr_coef);
    p_value_str_WM = num2str(p_value);
    
    
    h.Color = [1 1 1];
    h.Position = [1 1 1200 800];

    s = errorbar(median_median_R1_WM, median_median_diff_WM(1,:,pp), SE_diff_WM(1,:,pp), SE_diff_WM(1,:,pp), SE_R1_WM, SE_R1_WM,'o');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'r';
    s.MarkerFaceColor = [0.8 0.8 0.8];
    s.Color = 'r';
    hold on
    plot(median_median_R1_WM, median_median_diff_WM(1,:,pp),'ro');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_WM;
        y = median_median_diff_WM(1,:,pp);
        idx = isnan(y);
        y(idx) = [];
        x(idx) = [];
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
    ylabel(parameters{1,pp+1},'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{WM}=',corr_coef_str_WM,'; p=',p_value_str_WM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;

end
saveas(gcf,fullfile(outpath,plotname_WM),'fig');

% FDR correction using Benjamini-Hochberg method
if exist('mafdr', 'file') % Check if mafdr is available
    adjpvalues = mafdr(rawpvalues, 'BHFDR', true);
else
    % Manual FDR correction if mafdr is not available
    [sortedP, sortIdx] = sort(rawpvalues);
    numTests = length(pValues);
    adjustedP = sortedP .* numTests ./ (1:numTests)';
    adjustedP(adjustedP > 1) = 1; % Ensure p-values are not >1
    
    % Reorder back to original indices
    adjpvalues = zeros(size(rawpvalues));
    adjpvalues(sortIdx) = adjustedP;
end

length(correlations)
resultsTable = table(parameters(2:end)',correlations, rawpvalues, adjpvalues, ...
                     'VariableNames', {'Parameter','WM_Correlation', 'Raw_pValue', 'Adjusted_pValue'});
disp(resultsTable);

end