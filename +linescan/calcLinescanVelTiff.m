function Result = calcLinescanVelTiff(varargin)
    if nargin == 1 && istable(varargin{1})
        T = varargin{1};
    else
        p = parseInputs(varargin);
        T = getSettings(p);
    end
    
    %% Display settings table
    % TODO: always display all settings if any need to be set, but make
    % sure to set the default values to what was pulled from metadata
    % TODO: should input override metadata? Issue warning if don't match?
    T = uitableinteractive(T);
    
    %% Process Linescans
    % subtract average value of each column from image to take out vertical stripes
    % TODO: should this happen frame-by-frame, block-by-block or over
    % entire image??
    % Before this was happening block-by-block which seems a bit sus
    % Maybe should move out depening on which
    % TODO: Create DisplayWaitbar argument?
    % Create waitbar
    nFiles = size(T,1);
    
    % Create waitbar
    hWait = waitbar(0, '', 'Name', 'Linescan Files');
    hWait.Children.Title.Interpreter = 'none';
    for iFile = 1:1:nFiles
        % Update waitbar
        [~,name,~] = fileparts(T.filepath{iFile});
        waitbar((iFile-1)/nFiles, hWait, {['Processing linescan (', num2str(iFile), '/' num2str(nFiles),')'], name})
        
        filepath = T.filepath{iFile};
        MaxLines = T.MaxLines(iFile);
        UseAvg = T.UseAvg(iFile);
        % TODO: should these be rounded here?
        left = T.MaskLeft(iFile);
        right = T.MaskRight(iFile);
        msPerLine = T.msPerLine(iFile);
        umPerPx = T.umPerPx(iFile);
        WinSize = T.WinSize(iFile);
        WinStep = T.WinStep(iFile);
        transform = char(T.Transform{iFile});
        metric = char(T.Metric{iFile});
        optimizer = char(T.Optimizer{iFile});
        filtervar = T.FilterVar(iFile);
        useGPU = T.UseGPU(iFile);
        invert = T.Invert(iFile);
        
        % Read full image stack and reshape to 2D
        hTiffReader = util.io.readTiffStack(filepath);
        I = permute(hTiffReader.data(), [2,1,3]);
        delete(hTiffReader);

        nPix = size(I, 2);
        I = permute(I, [1,3,2]);
        I = reshape(I, [], nPix);
        nLines = size(I, 1);
        
        % TODO: validate this earlier?
        MaxLines = min(MaxLines, nLines);

        % Get rid of vertical stripes?
        if UseAvg
            I = double(I)-mean(I,1);
        end

        % Take requested section of image
        % TODO: need to round here? Should masklinescan return rounded
        % value?
        I = I(1:MaxLines, left:right);
        
        % Transfer image to GPU, if requested
        if useGPU; gpuDevice(1); I = gpuArray(I); end
        
        % Run Linescan
    %     tic
        Result = linescan.calcLinescanVel(I, msPerLine, umPerPx, WinSize, WinStep, 'Invert', invert, 'Transform', transform, 'Metric', metric, 'Optimizer', optimizer, 'FilterVar',filtervar);
    %     toc
        
        % TODO; allow user to specify this?
        Datafile = strrep(filepath,'.tif',[' rawVel ', num2str(WinStep), '-', num2str(WinSize), '-New.mat']);
        % TODO: make sure not missing any variables
        Settings = T(iFile, :);
        save(Datafile,'Settings','Result');
        
    end
    % Complete and close waitbar
    waitbar(1, hWait)
    close(hWait);
    
    if nargout < 1
        clear Result
    end
end

function p = parseInputs(args)
    p = inputParser();
    p.addOptional('filepath','',@ischar)
    p.addOptional('WinSize',75,@(x) isnumeric(x)&&isscalar(x));
    p.addOptional('WinStep',50,@(x) isnumeric(x)&&isscalar(x));
    % TODO: should these be parameters?
    p.addParameter('msPerLine',NaN,@(x) isnumeric(x)&&isscalar(x));
    p.addParameter('umPerPx',NaN,@(x) isnumeric(x)&&isscalar(x));
    p.addParameter('Mask','Visual',@(x) ischar(x)||isvector(x))
    p.addParameter('Transform','Radon',@ischar);
    p.addParameter('Metric','Var',@ischar);
    p.addParameter('Optimizer','radonlegacy',@ischar);
    p.addParameter('MaxLines',inf);
    p.addParameter('UseAvg',false,@islogical);
    p.addParameter('FilterVar',25,@isreal);
    p.addParameter('UseGPU',false,@islogical);
    p.addParameter('Invert',NaN,@islogical);
    p.parse(args{:});
end

function T = getSettings(p)
    if isempty(p.Results.filepath)
        % Get file to open
        [fname,pname] = uigetfile({'*.tif;*.tiff';'*.*'},'MultiSelect','on');
        Openfile = fullfile(pname, fname);
        % TODO: is this necessary?
        disp(Openfile);
    else
        Openfile = p.Results.filepath;
    end
    
    if ~iscell(Openfile)
        Openfile = {Openfile};
    end
    
    T = table();
    T.Properties.Description = 'Linescan Settings';
    for iFile = 1:1:length(Openfile)
        %% Metadata
        hTiff = Tiff(Openfile{iFile}, 'r');

        try software = hTiff.getTag('Software'); catch; software = ''; end
        imageDescription = hTiff.getTag('ImageDescription');
        try XResolution = hTiff.getTag('XResolution'); catch; XResolution = NaN; end
        try ResolutionUnit = hTiff.getTag('ResolutionUnit'); catch; ResolutionUnit = 1; end
        delete(hTiff);
        
        % Convert to microns
        switch ResolutionUnit
            case 1  % 'None'
                % TODO: is this a valid assumption?
                % Assumed to be um
                umPerPx = 1/XResolution;
            case 2  % 'Inch' = 25400 um
                umPerPx = 1/XResolution*25400;
            case 3  % 'Centimeter' = 10000 um
                umPerPx = 1/XResolution*10000;
        end
        if ~isnan(p.Results.umPerPx) && umPerPx ~= p.Results.umPerPx
            warning('Input umPerPx does not match XResolution in Tiff header data');
            % TODO: rename this variable
            umPerPx = p.Results.umPerPx;
        end

        % Determine which version of ScanImage from: 3.8, 2016b, 2018a
        % ScanImage 2016b, 2018a
        version = char(extractBetween(software,'SI.VERSION_MAJOR = ',newline));
        if isempty(version)
            % ScanImage 3.8
            version = char(extractBetween(imageDescription,'state.software.version=',char(13))); % char(13) is a newline character != newline
        end

        % Parse and format metadata for different ScanImage versions
        msPerLine = NaN;
        invert = false;
        switch version
            case '3.8'
                % ImageDescription tag stores 'state' structure
                s = util.io.str2struct(imageDescription);
                msPerLine = s.state.acq.msPerLine;
                invert = sign(s.state.acq.scanRotation) == -1;
            case {'2016','''2018a'''}
                s = util.io.str2struct(software);
                % TODO:
                % msPerLine = ?;
            case ''
                % Do nothing: no ScanImage version found
            otherwise
                % TODO: is this necessary? Set to NaN??
                warning(['Unsupported ScanImage version ', version, '. Cannot detect line period from metadata.'])
                % TODO: or display linked text to allow user to turn off this
                % warning
        end
        
        % Check if metadata matches input
        if ~isnan(p.Results.msPerLine) && msPerLine ~= p.Results.msPerLine
            warning('Input msPerLine does not match XResolution in Tiff header data');
            msPerLine = p.Results.msPerLine;
        end
        if ~isnan(p.Results.Invert) && invert ~= p.Results.Invert
            warning('Input invert does not match scanRotation in Tiff header data');
            invert = p.Results.Invert;
        end
        
        %% TODO: get mask
        % TODO: allow mask method to be set as input    % Mask linescan
        if ischar(p.Results.Mask)
            [MaskLeft, MaskRight] = linescan.maskLinescan(Openfile{iFile}, p.Results.Mask);
        else
            MaskLeft = p.Results.Mask(1);
            MaskRight = p.Results.Mask(2);
        end
        
        %% Display settings in table
        T.filepath{iFile} = Openfile{iFile};
        T.WinSize(iFile) = p.Results.WinSize;
        T.WinStep(iFile) = p.Results.WinStep;
        T.msPerLine(iFile) = msPerLine;
        T.umPerPx(iFile) = umPerPx;
        T.MaskLeft(iFile) = MaskLeft;
        T.MaskRight(iFile) = MaskRight;
        T.Transform{iFile} = categorical({p.Results.Transform}, {'Radon', 'Rotate'});
        T.Metric{iFile} = categorical({p.Results.Metric}, {'Var', 'Sep'});
        T.Optimizer{iFile} = categorical({p.Results.Optimizer}, {'globalsearch', 'multistart',...
        'binarysearch', 'exhaustive', 'radonlegacy', 'svdlegacy'});
        T.MaxLines(iFile) = p.Results.MaxLines;
        T.UseAvg(iFile) = p.Results.UseAvg;
        T.FilterVar(iFile) = p.Results.FilterVar;
        T.UseGPU(iFile) = p.Results.UseGPU;
        T.Invert(iFile) = invert;
    end
end