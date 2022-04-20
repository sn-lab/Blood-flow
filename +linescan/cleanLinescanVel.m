function Results = cleanLinescanVel(varargin)
%cleanLinescanVel Allows user to clean linescan velocity by setting upper
%and lower limits on velocity, then generates summary statistics.
%   Detailed explanation goes here

    p = inputParser();
    p.addOptional('filepath','',@ischar)
    p.parse(varargin{:});

    % TODO: make sure all output files are being saved in input directory;
    % Get input folder and file
    % TODO: allow multi-select
    if isempty(p.Results.filepath)
        % Get file to open
        [fname,pname] = uigetfile('*.*','MultiSelect','on');
        Openfile = fullfile(pname, fname);
    else
        Openfile = p.Results.filepath;
    end
    
    % For single file, put Openfile (char) into cell array to standardize
    % format
    if ~iscell(Openfile)
        Openfile = {Openfile};
    end
    
    % Intialize output table
    % TODO: preallocate size?
    Results = table();
    
    for iFile = 1:1:length(Openfile)
        % TODO: is this necessary?
        disp(Openfile);
        
        % TODO: deal with OpenFile being char vs cell array
        Data = load(Openfile{iFile});

        VelLimits = [-inf, inf];
        PctValid = 100;
        vel = Data.Result(:,3);
        avg = mean(Data.Result(:,3));
        stdev = std(Data.Result(:,3));
        
        time = Data.Result(:,2);
        
        % Create figure
        f = figure;
        set (f, 'Units', 'normalized', 'Position', [.1, .1, .8 .8], 'KeyPressFcn', @figureKeyPressFcn);
        
        % Plot velocity trace on left side of figure
        subplot('Position',[0.05, 0.1, 0.6, 0.8]);
        VelPlot = plot(time, vel, '-b.');
        title(Openfile{iFile}, 'Interpreter', 'none');
        subtitle(['Velocity: ', num2str(avg), 'mm/s  Std. Dev.: ', num2str(stdev)])
        xlabel('Time (s)')
        ylabel('Velocity (mm/s)')
        
        % Plot histogram on right side of figure
        subplot('Position', [0.7, 0.1, 0.25, 0.8]);
        % TODO: let MATLAB calculate number of bins?
        VelHist = histogram(vel,25,'Orientation','horizontal');
        
        % Wait for user to close the figure
        uiwait(f);
        
        % Save statistics in table
        T.filepath{iFile} = Openfile{iFile};
        T.Avg(iFile) = avg;
        T.Stdev(iFile) = stdev;
        T.PctValid(iFile) = PctValid;
        
        % Save cleaned data and settings, summary statistics in output file
        Data.Vel = vel;
        Data.VelLimits = VelLimits;
        Data.avg = avg;
        Data.stdev = stdev;
        Data.PctValid = PctValid;
        save(strrep(Openfile{iFile}, '.mat', '_ExcludedPts.mat'),'-struct','Data');
    end
    % Write CSV file with summary statistics for each input file
    filepath = fileparts(Openfile{1});
    writetable(Results, fullfile(filepath, 'Linescan Results.csv'));
    
    function figureKeyPressFcn(~,~)
        getNewVelLimits();
        updateData();
    end
    
    
    function updateData()
        % Deal with Outliers
        vel = Data.Result(:,3);
        outliers = (vel<=VelLimits(1)) | (vel>=VelLimits(2));
        vel(outliers) = NaN;
        
        % Recalculate summary statistics
        avg = mean(vel, 'omitnan');
        stdev = std(vel, 'omitnan');
        PctValid = 100*sum(~outliers)/length(vel);
        
        % Update plots
        VelPlot.YData = vel;
        VelHist.Data = vel;
        subtitle(['Velocity: ', num2str(avg), 'mm/s  Std. Dev.: ', num2str(stdev)])
    end
    
    function getNewVelLimits()
        % Prompt user for new upper and lower velocity limits
        limitok = false;
        while not(limitok) % check that upper limit > lower limit
            prompt = {'Upper Limit (mm/s)','Lower Limit (mm/s)'};
            definput = {num2str(VelLimits(2)), num2str(VelLimits(1))};
            dlgTitle='Velocity Limits';
            dims=1;
            opts = 'on';
            answer=inputdlg(prompt,dlgTitle,dims,definput,opts);
            UpLimit = str2double(answer(1));
            LowLimit = str2double(answer(2));

            limitok = UpLimit > LowLimit;
            VelLimits = [LowLimit, UpLimit];
        end %while not(limitok)
    end
    
end






            