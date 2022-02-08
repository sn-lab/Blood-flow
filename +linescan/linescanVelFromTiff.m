function Result = linescanVelFromTiff(varargin)
    p = inputParser();
    p.addOptional('filepath','',@ischar)
    p.addOptional('WinSize',75,@(x) isnumeric(x)&&isscalar(x));
    p.addOptional('WinStep',50,@(x) isnumeric(x)&&isscalar(x));
    % TODO: should these be parameters?
    p.addParameter('msPerLine',[],@(x) isnumeric(x)&&isscalar(x));
    p.addParameter('umPerPx',[],@(x) isnumeric(x)&&isscalar(x));
    p.addParameter('Mask','Visual',@(x) ischar(x)||isvector(x))
    p.addParameter('Method','Radon',@ischar);
    p.addParameter('Maxlines',inf);
    p.addParameter('UseAvg',false,@islogical);
    p.parse(varargin{:});
    
    
    % TODO: allow output as CSV?

    % TODO: make sure all output files are being saved in input directory;
    % Get input folder and file
    % TODO: allow multi-select
    if isempty(p.Results.filepath)
        % Get file to open
        [fname,pname] = uigetfile('*.*');
        Openfile = fullfile(pname, fname);
        disp(Openfile);
    else
        Openfile = p.Results.filepath;
    end
    
    %% Metadata
    hTiffReader = util.io.readTiffStack(Openfile);    
    meta = hTiffReader.metadata();
    imageDescription = hTiffReader.descriptions();
    delete(hTiffReader);
    
    hTiff = Tiff(Openfile, 'r');
    try
        model = hTiff.getTag('Model');
    catch
    end
    delete(hTiff);
    
    % TODO: get XResolution from Tiff tag instead?
    % TODO: need to put if model exists statement around this?
    s = util.io.str2struct(model);
    if s.MIPS.MIPS.Px2um.Active
        if s.MIPS.Px2um.VoxelSize(1) ~= s.MIPS.Px2um.VoxelSize(2)
            warning('Pixels are not square (different X and Y sizes). Using X size')
        end
        umPerPx = s.MIPS.Px2um.VoxelSize(1);
    end
        
    % Determine which version of ScanImage from: 3.8, 2016b, 2018a
    % TODO: Force user to use MIPS or change metadata in some way in order
    % to pull umPerPx from metadata. Too much to include full px2um package
    % here.

    %ScanImage 2016b, 2018a
    version = char(extractBetween(meta,'SI.VERSION_MAJOR = ',newline));
    if isempty(version)
        % ScanImage 3.8
        version = char(extractBetween(imageDescription{1,1},'state.software.version=',char(13))); % char(13) is a newline character != newline
    end
    
    % Parse and format metadata for different ScanImage versions
    switch version
        case '3.8'
            % 'state' structure stored in ImageDescription tag
            % No Software or Artist Tags
            s = util.io.str2struct(imageDescription{1,1});
            msPerLine = s.state.acq.msPerLine;
            
        case {'2016','''2018a'''}
            meta = strsplit(meta, [newline, newline]);
            software = meta{1};
            s = util.io.str2struct(software);
            % TODO:
            % msPerLine = ?;
        otherwise
            warning(['Unsupported ScanImage version ', version, '. Cannot detect line period from metadata.'])
    end
    % TODO: still need objective
    
    %% Settings
    % TODO: always display all settings if any need to be set, but make
    % sure to set the default values to what was pulled from metadata
    % TODO: should input override metadata? Issue warning if don't match?
    if isempty(p.Results.msPerLine)
        % Try to get msPerLine from tiff file, check if matches
        % If not, prompt user
    else
        msPerLine = p.Results.msPerLine;
    end
       
    if isempty(p.Results.umPerPx)
        % Try to get umPerPx from tiff file, check if matches
        % If not, prompt user
        % For now, build this from the list of versions in MIPS. In future,
        % might be good to use built-in ScanImage stuff
        
    else
        umPerPx = p.Results.umPerPx;
    end
    
    %%
    % subtract average value of each column from image to take out vertical stripes
    % TODO: should this happen frame-by-frame, block-by-block or over
    % entire image??
    % Before this was happening block-by-block which seems a bit sus
    % Maybe should move out depening on which
    
    % read in image data and reshape
    hTiffReader = util.io.readTiffStack(Openfile);
    I = permute(hTiffReader.data(), [2,1,3]);
    delete(hTiffReader);
    
    % Mask linescan
    if ischar(p.Results.Mask)
        [left, right] = linescan.maskLinescan(Openfile, p.Results.Mask);
    else
        left = p.Results.Mask(1);
        right = p.Results.Mask(2);
    end
    
    nPix = size(I, 2);
    I = permute(I, [1,3,2]);
    I = reshape(I, [], nPix);
    nLines = size(I, 1);
    
    Maxlines = min(p.Results.Maxlines, nLines);
    
    % Get rid of vertical stripes?
    if p.Results.UseAvg
        I = I-mean(I);
    end
    
    % Take requested section of image
    I = I(1:Maxlines,round(left):round(right));
    
    % Run Linescan
    errorcheck = false;
    WinSize = p.Results.WinSize;
    WinStep = p.Results.WinStep;
    method = p.Results.Method;
    Xfactor = umPerPx;
    Tfactor = 1/msPerLine;
    
%     tic
    Result = linescan.calcLinescanVel(I, msPerLine, umPerPx, WinSize, WinStep, errorcheck, 'Method', method);
%     toc
    
    % TODO: save Result
    % TODO: this should be moved out into calling function 
    Datafile = [char(strrep(Openfile,'.tif',[' rawVel ', num2str(WinStep), num2str(WinSize)])),'-New.mat'];
    save(Datafile,'Openfile','Result','WinSize','WinStep','Tfactor','Xfactor','left', 'right');

end