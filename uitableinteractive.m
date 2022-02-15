function T = uitableinteractive(T)
    % Create UI figure
    hUIFig = uifigure('Units','normalized','Position',[0.25, 0.25, 0.5, 0.5]);

    % Create table UI component
    hUITable = uitable(hUIFig, 'Data', T, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.9]);
    hUITable.ColumnSortable = true;
    hUITable.ColumnEditable = true;
    hUITable.DisplayDataChangedFcn = @updateData;
    hUITable.KeyPressFcn = @UITableKeyPressFcn;
    uiwait(hUIFig);

    function updateData(src,~)
        T = src.DisplayData;
    end
    

    function UITableKeyPressFcn(src,evt)
        % TODO: do this with a switch-case?
        if length(evt.Modifier) == 1 && strcmp(evt.Modifier, 'control') && strcmp(evt.Key, 'v')
            pasteTableData(src)
        elseif length(evt.Modifier) == 1 && strcmp(evt.Modifier, 'control') && strcmp(evt.Key, 'c')
            copyTableData(src)
        end
    end
    
    function pasteTableData(src)
        % TODO: what about pasting a scalar value into multiple selected
        % cells?
        % TODO: allow user to paste in more rows?
        
        data = src.DisplayData;
        str = clipboard('paste');

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
            % TODO: need to handle other data classes?
            % Need to access in this way to get true type, otherwise will
            % return as cell
            % TODO: maybe return everything as a acell array and then index
            % the cell array to get the class...
            % dataCol = data.(data.Properties.VariableNames{i});
            dataCol = data{:,i};
            if iscell(dataCol); dataCol = dataCol{:}; end
            switch class(dataCol)
                case 'char'
                    % TODO: really this should use curly braces
                    data{toRows,toCols(i)} = cellstr(fromRows,fromCols(i));
                case 'cell'
                    data(toRows,toCols(i)) = cellstr(fromRows,fromCols(i));
                case 'double'
                    data{toRows,toCols(i)} = str2double(cellstr(fromRows,fromCols(i)));
                case 'logical'
                    data{toRows,toCols(i)} = cellfun(@(x) strcmpi(x,'true'), cellstr(fromRows,fromCols(i)));
                case 'categorical'
                    data{toRows,toCols(i)} = {categorical(cellstr(fromRows,fromCols(i)),categories(dataCol))};
                otherwise
                    error(['Cannot paste data into column ', num2str(i),...
                           'due to unsupported type: ', class(dataCol)])
            end
        end
        % TODO: why setting src.Data rather than src.DisplayData?
        src.Data = data;
        updateData(src)
    end

    function copyTableData(src)
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
        % TODO: is it necessary to strip trailing newline?
        str = str(1:end-1);
        
        % Set computer clipboard
        clipboard('copy', str);
    end

end