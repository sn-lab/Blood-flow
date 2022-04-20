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
    
    vars = repmat({pname,'',75,50,NaN,NaN,0,1024,'Radon',inf,false}, nFiles, 1);
    vars{:,2} = fname;
    vars = [vars;vars;vars;vars];
    
    T = cell2table(vars,'VariableNames',{'Folder','File','WinSize','WinStep','msPerLine','umPerPx','MaskLeft','MaskRight','Method','Maxlines','UseAvg'});
    
%     assignin('base', 'T', T);
%     openvar('T')
%     disp('Check settings and press any key to continue')
%     pause
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
            pasteTableData(src)
        elseif length(evt.Modifier) == 1 && strcmp(evt.Modifier, 'control') && strcmp(evt.Key, 'c')
            copyTableData(src)
        end
    end
    
    function pasteTableData(src)
        % TODO: move this conent to separate function for better
        % organization?
        data = src.DisplayData;
        str = clipboard('paste');
        % TODO: maybe could combine more of the if/else clauses
        % Should this go row-by-row, column-by-column or cell by cell?
        % Column-by-column seems to make most sense because datatype is
        % consistent

        % Parse clipboard string
        str = strtrim(str);
        cellstr = splitlines(str);
        cellstr = cellfun(@(x) strsplit(x, '\t'), cellstr, 'UniformOutput', false);
        cellstr = vertcat(cellstr{:});

        % Indexes
        % TODO: consolidate this once it's tested
        toRowStart = src.Selection(1,1);
        toRowEnd = min(src.Selection(1)+size(cellstr,1)-1, size(data,1));
        toRows = toRowStart:toRowEnd;

        fromRowStart = 1;
        fromRowEnd = 1+length(toRows)-1;
        fromRows = fromRowStart:fromRowEnd;

        toColStart = src.Selection(1,2);
        toColEnd = min(size(cellstr, 2), size(data, 2));
        toCols = toColStart:toColEnd;

        fromColStart = 1;
        fromColEnd = 1+length(toCols)-1;
        fromCols = fromColStart:fromColEnd;

        for i = 1:1:length(toCols)
            switch class(data{:,toCols(i)})
                case 'cell'
                    data(toRows,toCols(i)) = cellstr(fromRows,fromCols(i));
                case 'double'
                    data{toRows,toCols(i)} = str2double(cellstr(fromRows,fromCols(i)));
                case 'logical'
                    data{toRows,toCols(i)} = cellfun(@(x) strcmpi(x,'true'), cellstr(fromRows,fromCols(i)));
                otherwise
                    error(['Cannot paste data into column ', num2str(col),...
                           'due to unsupported type: ', class(data{:,col})])
            end
        end
        src.Data = data;
    end

    function copyTableData(src)
        % Add newlines
        % Run cellfun to add tabs?
        % Probably need for-loop from min(row) to max(row) and min(col) to
        % max(col) in order to add apropriate number of tabs etc. for
        % non-rectangular selections
        % TODO: force rectangular selections??
        % TODO: maybe faster to take rectangular selection and then replace
        % non-selected cells with empty strings, then concatenate
        % everything with tabs etc
        % First, concatenate rows together to make vector, then do rows
        cellstr = table2cell(src.DisplayData);
        cellstr = string(cellstr);
        
        % Set all unselected cells to empty
        selected = src.Selection;
        sz = size(src.DisplayData);
        selectedInd = sub2ind(size(src.DisplayData), selected(:,1), selected(:,2));
        unselectedTF = true(sz);
        unselectedTF(selectedInd) = false;
        cellstr(unselectedTF) = '';
        
        % Take rectangle around selected cells, with unselected cells set
        % to empty
        starts = min(selected);
        ends = max(selected);
        cellstr = cellstr(starts(1):ends(1), starts(2):ends(2));
        
        % Format as tab/newline delimited string
        str = '';
        for iRow = 1:1:size(cellstr, 1)
            rowStr = sprintf('%s\t',cellstr(iRow,:));
            str = sprintf('%s%s\n',str, rowStr(1:end-1));
        end
        % TODO: is this necessary?
        str = str(1:end-1);
        
        % Now take subset of min to max
        
        clipboard('copy', str);
    end

end