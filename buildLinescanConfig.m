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
    
    
%     for iFile = 1:1:nFiles
%         % Mask linescan
%         Openfile = fullfile(T.Folder{iFile}, T.File{iFile});
%         [T.MaskLeft(iFile), T.MaskRight(iFile)] = linescan.maskLinescan(Openfile, 'Visual');
%     end
    
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
    hUITable.KeyPressFcn = @UITableKeyPressFcn;
    uiwait(hUIFig);

        function updateData(~,~)
            T = hUITable.DisplayData;
        end
    
    function UITableKeyPressFcn(src,evt)
%         data = get(handles.uitable1,'Data');
        if length(evt.Modifier) == 1 && strcmp(evt.Modifier, 'control') && strcmp(evt.Key, 'v')
            data = src.DisplayData;
            clipboard = importdata('-pastespecial');
            clipboard = cellfun(@(c) strsplit(c,'\t','CollapseDelimiters',false), clipboard, 'UniformOutput', false);
            for i = 1:1:length(clipboard)
                row = i+src.Selection(1)-1;
                cols = src.Selection(1):src.Selection(1)+length(clipboard{i,1})-1;
                data(row,cols) = clipboard{i,1};
            end
            src.DisplayData = data;
        end

    end
end