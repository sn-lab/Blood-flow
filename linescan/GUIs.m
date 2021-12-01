    % Calibration factors
    % Calculate Tfactor (number of pixels per ms)
    % TODO: try to get this info from the image metadata
    prompt={'ms per line', 'microns per pixel'};
    def={'0.6', '0.1'};
    dlgTitle='Conversion factors';
    lineNo=1;
    answer2=inputdlg(prompt,dlgTitle,lineNo,def,'on');
    Tfactor = 1/str2double(cell2mat(answer2(1))); % ypixel per ms
    Xfactor = str2double(cell2mat(answer2(2))); % microns per xpixel
    
    % Running parameters from user
prompt={'Number of time pixels per data point', 'time pixels between data points',...
    'number of lines to process', 'Start with file #:'};
def={'75', '50', '50000', '1', 'no'};
dlgTitle='Processing parameters';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def,'on');
WinSize =  str2double(answer(1));   % actual data used is only center circle ~70% of area (square window)
WinPixelsDown = str2double(answer(2)); % number of pixels between top of last window and next window
Maxlines =  str2double(answer(3)); %  total number of lines
startfilenumber = str2double(answer(4));


answer = questdlg('Turn on display?', 'Debugging Mode', 'Yes', 'No', 'No');
if strcmp(answer, 'No')
    errorcheck = 0;
else
    errorcheck = 1;
end