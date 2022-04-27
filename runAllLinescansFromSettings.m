%% 
folder = '/home/schafferlab/Desktop/Nancy/LongLinescan/All Data/Images/CH4';

T = readtable('/home/schafferlab/Desktop/Nancy/LongLinescan/LinescanSettings.xlsx');
T.filepath = fullfile(folder,T.filepath);
T.MaxLines = str2double(T.MaxLines);
T = T(99:end, :);
linescan.calcLinescanVelTiff(T);

%%
T = readtable('/mnt/SN Data Server 2/Nancy Ruiz/LongLinescan/LinescanSettings.xlsx');
folder = '/mnt/SN Data Server 2/Nancy Ruiz/LongLinescan/All Data/Images/CH4';
T.filepath = fullfile(folder,T.filepath);
rows = [130, 163, 164, 165, 167, 168, 242, 243];
T = T(rows,:);
T.MaxLines = str2double(T.MaxLines);
linescan.calcLinescanVelTiff(T);