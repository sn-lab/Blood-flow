function linescanVelFromTiff(varargin)
    p = inputParser();
    p.addOptional(filepath,'',@ischar)
    p.addParameter('msPerLine',[],@(x) isnumeric(x)&&isscalar(x));
    p.addParameter('umPerPx',[],@(x) isnumeric(x)&&isscalar(x));
    p.addOptional('Mask','Auto',@ischar)
     % TODO: these should probably be moved up to calling function
    p.addParameter('Method','Radon',@ischar);
    p.addParameter('Maxlines',inf);
    p.addParameter('UseAvg',false,@islogical);
    p.parse(varargin{:});
    
    
    % TODO: allow output as CSV?

    % TODO: make sure all output files are being saved in input directory;
    % Get input folder and file
    % TODO: allow multi-select
    if isempty(p.Results.Openfile)
        % Get file to open
        [fname,pname] = uigetfile('*.*');
        Openfile = fullfile(pname, fname);
        disp(Openfile);
    else
        Openfile = p.Results.Openfile;
        [pname,fname] = fileparts(Openfile);
    end
    
    
    if isempty(p.Results.msPerLine)
        % Try to get msPerLine from tiff file
        % If not, prompt user
    end
    
    if isempty(p.Results.umPerPx)
        % Try to get umPerPx from tiff file
        % If not, prompt user
    end
    
    Maxlines = min(Maxlines, nLines);
    
    % subtract average value of each column from image to take out vertical stripes
    % TODO: should this happen frame-by-frame, block-by-block or over
    % entire image??
    % Before this was happening block-by-block which seems a bit sus
    % Maybe should move out depening on which
    if UseAvg
        I = I-mean(I);
    end
    
     % read in image data and reshape
    I = f_loadTiff(Openfile);
    nPix = size(I, 2);
    I = permute(I, [1,3,2]);
    I = reshape(I, [], nPix);
    nLines = size(I, 1);
    
    if ischar(p.Results.Mask)
        [left, right] = maskLinescan(I,p.Results.Mask);
    else
        left = p.Results.Mask(1);
        right = p.Results.Mask(2);
    end
    % TODO: apply mask
    
end


% TODO: probably just use tiff reader instead
function I = f_loadTiff(filename)
    T = Tiff(filename, 'r');
    
    % Extracting Matlab Type Name and Bits per Sample
    typeName = T.getTag('SampleFormat');
    if typeName == 1
       typeName = 'uint';
    elseif typeName  == 2
       typeName = 'int';
    else
       error('Tiff Sample Format is not supported. It must be UInt or Int');
    end
    bits = T.getTag('BitsPerSample');
    type = [typeName, num2str(bits)];

    % Get Image Height and Width
    width = T.getTag('ImageWidth');
    height = T.getTag('ImageLength');
            
    % Get number of images in Tiff Stack
    T.setDirectory(1);
    numImages = 1;
    while ~T.lastDirectory
       numImages = numImages + 1;
       T.nextDirectory();
    end

    I = zeros(height, width, numImages, type);
    % Setting data for each Image File Directory
    for IFD = 1:1:numImages
       T.setDirectory(IFD);
       I(:,:,IFD) = T.read();
    end
        
    close(T);
end