% For masking several files and storing values 'left' and 'right'
folder = '/Volumes/Extreme SSD/LongLinescan/All Data/Images/CH4/';

%% 
files = dir(fullfile(folder, '*.tif'));
diary('Mask Log.txt')
for i = 1:1:length(files)
    filename = fullfile(folder, files(i).name)
    [left, right] = linescan.maskLinescan(filename, 'Visual')
    files(i).Left = left;
    files(i).Right = right;
  end

%%
save('Masks.mat', 'files');