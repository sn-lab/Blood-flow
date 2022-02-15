function Result = calcLinescanVelTiff(varargin)
    p = inputParser();
    p.addOptional('filepath','',@ischar)
    p.addOptional('WinSize',75,@(x) isnumeric(x)&&isscalar(x));
    p.addOptional('WinStep',50,@(x) isnumeric(x)&&isscalar(x));
    % TODO: should these be parameters?
    p.addParameter('msPerLine',NaN,@(x) isnumeric(x)&&isscalar(x));
    p.addParameter('umPerPx',NaN,@(x) isnumeric(x)&&isscalar(x));
    p.addParameter('Mask','Visual',@(x) ischar(x)||isvector(x))
    p.addParameter('Method','Radon',@ischar);
    p.addParameter('Optimizer','fminbnd',@ischar);
    p.addParameter('MaxLines',inf);
    p.addParameter('UseAvg',false,@islogical);
    p.parse(varargin{:});
    
    
    % TODO: allow output as CSV?

    % TODO: make sure all output files are being saved in input directory;
    % Get input folder and file
    % TODO: allow multi-select
    if isempty(p.Results.filepath)
        % Get file to open
        [fname,pname] = uigetfile('*.*','MultiSelect','on');
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
                umPerPx = NaN;
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
        switch version
            case '3.8'
                % ImageDescription tag stores 'state' structure
                s = util.io.str2struct(imageDescription);
                msPerLine = s.state.acq.msPerLine;
            case {'2016','''2018a'''}
                s = util.io.str2struct(software);
                % TODO:
                % msPerLine = ?;
            otherwise
                % TODO: is this necessary? Set to NaN??
                warning(['Unsupported ScanImage version ', version, '. Cannot detect line period from metadata.'])
                % TODO: or display linked text to allow user to turn off this
                % warning
        end
        
        % Check if metadata matches input
        if ~isnan(p.Results.msPerLine) && msPerLine ~= p.Results.msPerLine
            warning('Input msPerLine does not match XResolution in Tiff header data');
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
        T.Method{iFile} = categorical({p.Results.Method}, {'Radon', 'SVD'});
        T.Optimizer{iFile} = categorical({p.Results.Optimizer}, {'fminbnd', 'legacy'});
        T.MaxLines(iFile) = p.Results.MaxLines;
        T.UseAvg(iFile) = p.Results.UseAvg;
        
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
    for iFile = 1:1:length(Openfile)
        MaxLines = T.MaxLines(iFile);
        UseAvg = T.UseAvg(iFile);
        % TODO: should these be rounded here?
        left = T.MaskLeft(iFile);
        right = T.MaskRight(iFile);
        msPerLine = T.msPerLine(iFile);
        umPerPx = T.umPerPx(iFile);
        WinSize = T.WinSize(iFile);
        WinStep = T.WinStep(iFile);
        errorcheck = false;
        method = char(T.Method{iFile});
        optimizer = char(T.Optimizer{iFile});
        
        % Read full image stack and reshape to 2D
        hTiffReader = util.io.readTiffStack(Openfile{iFile});
        I = permute(hTiffReader.data(), [2,1,3]);
        delete(hTiffReader);

        nPix = size(I, 2);
        I = permute(I, [1,3,2]);
        I = reshape(I, [], nPix);
        nLines = size(I, 1);

        MaxLines = min(MaxLines, nLines);

        % Get rid of vertical stripes?
        if UseAvg
            I = I-mean(I);
        end

        % Take requested section of image
        % TODO: need to round here? Should masklinescan return rounded
        % value?
        I = I(1:MaxLines, left:right);

        % Run Linescan
    %     tic
        Result = linescan.calcLinescanVel(I, msPerLine, umPerPx, WinSize, WinStep, errorcheck, 'Method', method, 'Optimizer', optimizer);
    %     toc

        % TODO: save Result
        % TODO: this should be moved out into calling function 
        Datafile = [char(strrep(Openfile,'.tif',[' rawVel ', num2str(WinStep), num2str(WinSize)])),'-New.mat'];
        % TODO: make sure not missing any variables
        save(Datafile,'Openfile','Result','WinSize','WinStep','msPerLine','umPerPx','MaxLines','UseAvg','left','right','method','optimizer');
    end
end