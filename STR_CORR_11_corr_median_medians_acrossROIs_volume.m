function [corr_acrossROIs_volume_GM,corr_acrossROIs_volume_WM,rho_pcorr,pval_pcorr] = STR_CORR_11_corr_median_medians_acrossROIs_volume(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,GM_ROI_vol,WM_ROI_vol,GM_thickness,GM_MeanR1,GM_GrayVol,GM_MeanCurv,outpath,confidence_bounds)

% calculates correlation coefficient matrix computing the median of the
% medians across subjects and creates scatterplots, after removing outliers
% for each parameter
N_par = length(parameters);
N_GM_ROIs = size(median_R1_GM,2);
N_WM_ROIs = size(median_R1_WM,2);
N_subjs = size(median_R1_GM,1);

corr_acrossROIs_volume_GM = zeros(N_par-1,2);
corr_acrossROIs_volume_WM = zeros(N_par-1,2);

plotname_GM = 'Scatterplot_GM_ROIvolume_vs_R1_acrossROIs';
plotname_WM = 'Scatterplot_WM_ROIvolume_vs_R1_acrossROIs';

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
    out_idx = isoutlier(tmp); % outliers = above/below 3 scaled deviations = median absolute deviations
    median_median_diff_GM(1,out_idx,pp) = NaN;
end

for pp = 1 : N_par-1
    tmp = median_median_diff_WM(1,:, pp);
    out_idx = isoutlier(tmp);
    median_median_diff_WM(1,out_idx,pp) = NaN;
end

SE_R1_GM = mad_median_R1_GM.*1.4826./sqrt(N_subjs);
SE_diff_GM = mad_median_diff_GM.*1.4826./sqrt(N_subjs);

SE_R1_WM = mad_median_R1_WM.*1.4826./sqrt(N_subjs);
SE_diff_WM = mad_median_diff_WM.*1.4826./sqrt(N_subjs);

correlations = zeros(N_par-1,1);
rawpvalues = zeros(N_par-1,1);

PositionArray = [1 1 350 400];
%% ROI volume
h = figure;

% median across subjects for each region
       
[r,p] = corrcoef(median_median_R1_GM, GM_ROI_vol(:,1), 'rows','complete');
size(median_median_R1_GM)

correlations(1,1) = r(2);
rawpvalues(1,1) = p(2);

    corr_acrossROIs_volume_GM(1,1) = r(2);
    corr_acrossROIs_volume_GM(1,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];
    h.Position = PositionArray;


    plot(median_median_R1_GM, GM_ROI_vol(:,1),'bo');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_GM;
        y = GM_ROI_vol(:,1);
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
    ylabel("GM ROI volume (Nvoxels)",'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;


saveas(gcf,fullfile(outpath,plotname_GM),'fig');
%% Thickness
h = figure;

% median across subjects for each region
       
[r,p] = corrcoef(median_median_R1_GM, GM_thickness, 'rows','complete');
size(median_median_R1_GM)

correlations(1,1) = r(2);
rawpvalues(1,1) = p(2);

    corr_acrossROIs_volume_GM(1,1) = r(2);
    corr_acrossROIs_volume_GM(1,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];
    h.Position = PositionArray;


    plot(median_median_R1_GM, GM_thickness,'bo');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_GM;
        y = GM_thickness;
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
    ylabel("GM thickness (mm)",'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;

plotname_GM = 'Scatterplot_GM_Thickness_vs_R1_acrossROIs';
saveas(gcf,fullfile(outpath,plotname_GM),'fig');
%{
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
%}
%% MeanCurv
h = figure;

% median across subjects for each region
       
[r,p] = corrcoef(median_median_R1_GM, GM_MeanCurv, 'rows','complete');
size(median_median_R1_GM)

correlations(1,1) = r(2);
rawpvalues(1,1) = p(2);

    corr_acrossROIs_volume_GM(1,1) = r(2);
    corr_acrossROIs_volume_GM(1,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];
    h.Position = PositionArray;


    plot(median_median_R1_GM, GM_MeanCurv,'bo');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_GM;
        y = GM_MeanCurv;
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
    ylabel("GM curvature (mm^{-1})",'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;

plotname_GM = 'Scatterplot_GM_Curvature_vs_R1_acrossROIs';
saveas(gcf,fullfile(outpath,plotname_GM),'fig');

%% GM R1
h = figure;

% median across subjects for each region
       
[r,p] = corrcoef(median_median_R1_GM, GM_MeanR1, 'rows','complete');
size(median_median_R1_GM)

correlations(1,1) = r(2);
rawpvalues(1,1) = p(2);

    corr_acrossROIs_volume_GM(1,1) = r(2);
    corr_acrossROIs_volume_GM(1,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];
    h.Position = PositionArray;


    plot(median_median_R1_GM, GM_MeanR1,'bo');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_GM;
        y = GM_MeanR1;
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
    ylabel("GM R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;

plotname_GM = 'Scatterplot_GM_R1_vs_R1_acrossROIs';
saveas(gcf,fullfile(outpath,plotname_GM),'fig');
%% GM volume by Freesurfer
h = figure;

% median across subjects for each region
       
[r,p] = corrcoef(median_median_R1_GM, GM_GrayVol, 'rows','complete');
size(median_median_R1_GM)

correlations(1,1) = r(2);
rawpvalues(1,1) = p(2);

    corr_acrossROIs_volume_GM(1,1) = r(2);
    corr_acrossROIs_volume_GM(1,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];
    h.Position = PositionArray;


    plot(median_median_R1_GM, GM_GrayVol,'bo');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_GM;
        y = GM_GrayVol;
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
    ylabel("GM GrayVol (mm^{3})",'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;

plotname_GM = 'Scatterplot_GM_Vol_vs_R1_acrossROIs';
saveas(gcf,fullfile(outpath,plotname_GM),'fig');
%% WHITE MATTER


correlations = zeros(N_par-1,1);
rawpvalues = zeros(N_par-1,1);

h = figure;

% median across subjects for each region
       
[r,p] = corrcoef(median_median_R1_WM, WM_ROI_vol(:,1), 'rows','complete');

    correlations(1,1) = r(2);
    rawpvalues(1,1) = p(2);

    corr_acrossROIs_noout_WM(1,1) = r(2);
    corr_acrossROIs_noout_WM(1,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_WM = num2str(corr_coef);
    p_value_str_WM = num2str(p_value);
    
    
    h.Color = [1 1 1];
    h.Position = [1 1 650 600] ;

 
    plot(median_median_R1_WM, WM_ROI_vol(:,1),'ro');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_R1_WM;
        y = WM_ROI_vol(:,1);
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
    ylabel("WM ROI volume (Nvoxels)",'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{WM}=',corr_coef_str_WM,'; p=',p_value_str_WM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;


saveas(gcf,fullfile(outpath,plotname_WM),'fig');

%% correlation bw fneurite and WM ROI vol
plotname_WM = 'Scatterplot_WM_ROIvolume_vs_fneurite_acrossROIs';
h = figure;

% median across subjects for each region
       size(median_median_diff_WM)
       size(WM_ROI_vol)

[r,p] = corrcoef(median_median_diff_WM(1,:,2), WM_ROI_vol(:,3), 'rows','complete');

    correlations(1,1) = r(2);
    rawpvalues(1,1) = p(2);

    corr_acrossROIs_noout_WM(1,1) = r(2);
    corr_acrossROIs_noout_WM(1,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_WM = num2str(corr_coef);
    p_value_str_WM = num2str(p_value);
    
    
    h.Color = [1 1 1];
    h.Position = [1 1 650 600];

 
    plot(median_median_diff_WM(1,:,2), WM_ROI_vol(:,3),'ro');


    if confidence_bounds == 1
        % confidence bounds
        x = median_median_diff_WM(1,:,2);
        y = WM_ROI_vol(:,3);
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
    ylabel("WM ROI volume (Nvoxels)",'FontSize',14,'FontWeight','Bold');
    xlabel("fneurite",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{WM}=',corr_coef_str_WM,'; p=',p_value_str_WM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;


saveas(gcf,fullfile(outpath,plotname_WM),'fig');

%% evaluates partial correlations in WM 
disp('fneurite vs R1 in WM controlling for Nvoxels')
x = squeeze(median_median_diff_WM(1,:,2));
x = x';
y = median_median_R1_WM';
z = WM_ROI_vol(:,1);
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

disp('Din vs R1 in WM controlling for Nvoxels')
x = squeeze(median_median_diff_WM(1,:,6));
x = x';
y = median_median_R1_WM';
z = WM_ROI_vol(:,1);
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

%% evaluates partial correlations in GM 
disp('fneurite vs R1 in GM controlling for Nvoxels, cortical thickness and curvature')
x = squeeze(median_median_diff_GM(1,:,2));
x = x';
y = median_median_R1_GM';
z = [GM_ROI_vol(:,1) GM_thickness GM_MeanCurv];
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

disp('fsoma vs R1 in GM controlling for Nvoxels, cortical thickness and curvature')
x = squeeze(median_median_diff_GM(1,:,3));
x = x';
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

disp('fextra vs R1 in GM controlling for Nvoxels, cortical thickness and curvature')
x = squeeze(median_median_diff_GM(1,:,4));
x = x';
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

disp('Rsoma vs R1 in GM controlling for Nvoxels, cortical thickness and curvature')
x = squeeze(median_median_diff_GM(1,:,5));
x = x';
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

disp('Din vs R1 in GM controlling for Nvoxels, cortical thickness and curvature')
x = squeeze(median_median_diff_GM(1,:,6));
x = x';
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

disp('De vs R1 in GM controlling for Nvoxels, cortical thickness and curvature')
x = squeeze(median_median_diff_GM(1,:,7));
x = x';
[rho_pcorr, pval_pcorr] = partialcorr(x,y,z,'rows','complete');
% Display results
fprintf('Partial correlation: %.4f\n', rho_pcorr);
fprintf('p-value: %.4f\n', pval_pcorr);

end