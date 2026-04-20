function [V_R1_maps, V_diff_maps] = STR_CORR_03_getting_metrics(initpath,parameters_list)
%STR_CORR_03_getting_metrics gets the R1 and diffusion metrics from the
%given paths
% Directory with subject data
cd(initpath)
folders = dir('pil*');
N_subj = length(folders);  % Number of subjects
N_par = length(parameters_list); % Number of diffusion parameters
% cell initialization
V_R1_maps = cell(N_subj,1);
V_diff_maps = cell(N_subj,N_par-1);

% getting R1, FR, SANDI
for ii = 1 : N_subj  
    
    disp('Filling up cells with R1 and diffusion metrics')
    disp(ii)

    for pp = 1 : length(parameters_list)

        par = parameters_list{1,pp}; 

%         if strcmp(par,'FR')
%             img_path_SANDI = strcat('/storage/ekaterina/Charmed/Dati_newFR/',folders(ii).name,'/',par,'2std.nii.gz');
%             nii = load_untouch_nii(img_path_SANDI);
%             V_SANDI_tot = double(nii.img);
%             V_diff_maps{ii,1} = V_SANDI_tot;
        if strcmp(par,'fw')
            img_path_SANDI = strcat('/storage/ekaterina/Charmed/Dati_SANDI/',folders(ii).name,'/',par,'2std.nii.gz');
            nii = load_untouch_nii(img_path_SANDI);
            V_SANDI_tot = double(nii.img);
            V_diff_maps{ii,1} = V_SANDI_tot;
        elseif strcmp(par,'R1')
            img_path_R1 = strcat('/storage/ekaterina/Charmed/Dati_R1/',folders(ii).name,'/',par,'2std.nii.gz');
            nii = load_untouch_nii(img_path_R1);
            V_R1_tot = double(nii.img);
            V_R1_maps{ii,1} = V_R1_tot;
        else
            img_path_SANDI = strcat('/storage/ekaterina/Charmed/Dati_SANDI/',folders(ii).name,'/',par,'2std.nii.gz');
            nii = load_untouch_nii(img_path_SANDI);
            V_SANDI_tot = double(nii.img);
            V_diff_maps{ii,pp-1} = V_SANDI_tot;
        end
   
    end
end
    % Final Debug Check
    disp(['Final V_R1_maps type: ', class(V_R1_maps)]);
    if ~isempty(V_R1_maps)
        disp(['First element type: ', class(V_R1_maps{1})]);
    end
end