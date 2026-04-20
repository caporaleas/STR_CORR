function done = save_subfolders_list_to_txt(mainFolder,outpath)
%save_subfolders_list_to_txt creates a txt file with the list of subfolders

% Get a list of all items in the folder
cd(mainFolder)
items = dir('pil*');

% Filter only directories whose names start with "pil"
subfolders = {items([items.isdir]).name};
subfolders = subfolders(startsWith(subfolders, 'pil'));

% Define the output text file path
outputFile = fullfile(outpath, 'subjects_list.txt');

% Write the list of subfolders to a text file
fileID = fopen(outputFile, 'w');
if fileID == -1
    error('Could not open file for writing.');
end

fprintf(fileID, '%s\n', subfolders{:});
fclose(fileID);

done = 'List of subjects saved to subjects_list.txt';
disp(done);

cd(outpath)

end