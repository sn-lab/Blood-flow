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
    
    % subtract average value of each column from image to take out vertical stripes
    % TODO: should this happen frame-by-frame, block-by-block or over
    % entire image??
    % Before this was happening block-by-block which seems a bit sus
    % Maybe should move out depening on which
    
    % read in image data and reshape
    h = util.io.readTiffStack(Openfile);
    I = permute(h.data(), [2,1,3]);
    delete(h);
    
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