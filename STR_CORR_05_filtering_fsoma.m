function [V_R1_maps_filtered,V_diff_maps_filtered] = STR_CORR_05_filtering_fsoma(initpath,fsomafilter,fsoma_min,V_R1_maps,V_diff_maps,parameters)
%APPLIES fsoma filtering both in R1 and diffusion maps (eliminates voxels
%where fsoma is lower than the minimum 
cd(initpath)
folders = dir('pil*');
N_subj = 20;  % Number of healthy subjects

if fsomafilter == 1
    for ii = 1 : N_subj  
        fsoma = V_diff_maps{ii,3};
        R1tmp = V_R1_maps{ii,1};
        R1tmp(fsoma < fsoma_min) = 0;
        V_R1_maps{ii,1} = R1tmp;
        for pp = 1 : length(parameters)-1
            tmp = V_diff_maps{ii,pp};
            tmp(fsoma < fsoma_min) = 0;
            V_diff_maps{ii,pp} = tmp;
        end
       
    end
    V_R1_maps_filtered = V_R1_maps;
    V_diff_maps_filtered = V_diff_maps;
    disp('fsoma filter applied')
else
    V_R1_maps_filtered = V_R1_maps;
    V_diff_maps_filtered = V_diff_maps;
    disp('No filter on fsoma applied')
end
end