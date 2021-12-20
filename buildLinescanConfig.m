function T = buildLinescanConfig()
    % TODO: allow output as CSV?

    % TODO: make sure all output files are being saved in input directory;
    % Get input folder and file
    % TODO: allow multi-select

    [fname,pname] = uigetfile('*.*','MultiSelect','on');
    
    
    if iscell(fname)
        nFiles = length(fname);
    else
        nFiles = 1;
    end
    
    vars = repmat({pname,'',75,50,[],[],0,1024,'Radon',inf,false}, nFiles, 1);
    vars(:,2) = fname;
    
    T = cell2table(vars,'VariableNames',{'Folder','File','WinSize','WinStep','msPerLine','umPerPx','MaskLeft','MaskRight','Method','Maxlines','UseAvg'});
    
    
    for iFile = 1:1:nFiles
        % Mask linescan
        Openfile = fullfile(T.Folder{iFile}, T.File{iFile});
        [T.MaskLeft(iFile), T.MaskRight(iFile)] = linescan.maskLinescan(Openfile, 'Visual');
    end
    
    T = displaySettings(T);
end


function T = displaySettings(T)
    % Create table array

    % Create UI figure
    hUIFig = uifigure('Units','normalized','Position',[0.25, 0.25, 0.5, 0.5]);

    % Create table UI component
    hUITable = uitable(hUIFig, 'Data', T, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.9]);
    hUITable.ColumnSortable = true;
    hUITable.ColumnEditable = true;
    hUITable.DisplayDataChangedFcn = @updateData;
    uiwait(hUIFig);

        function updateData(~,~)
            T = hUITable.DisplayData;
        end
end