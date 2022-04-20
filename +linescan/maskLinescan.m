function [WinLeft, WinRight] = maskLinescan(Openfile, method)
    switch method
        case 'Visual'
            [WinLeft, WinRight] = maskLinescanVisual(Openfile);
        case 'Auto'
            [WinLeft, WinRight] = maskLinescanAuto(I);
            error('Not supported');
        otherwise
            error('Unrecognized method');
    end
end

function [WinLeft, WinRight] = maskLinescanAuto(I)
    proj = max(I, [], [2,3]);
    proj = movmean(proj, round(length(proj)/32));
    % TODO: allow user to pass pctMax?
    TF = proj > 0.4*max(proj);
    dTF = diff(TF);
    rise = find(dTF == 1);
    fall = find(dTF == -1);
    
    % Get number of widest peak
    w = fall-rise;
    [m,iMaxW] = max(w);
    
    if m < 100
        warning('Detected window width is small (<100 px). Detection may be inaccurate')
    end
    
    WinLeft = rise(iMaxW);
    WinRight = fall(iMaxW);
    %TODO: add an option for displaying the detected region
%     figure; imshow(I(:,:,1), []); drawrectangle('Position', [WinLeft, 1, WinRight-WinLeft, size(I,1)],'Color','r');
end
    

    
function [WinLeft, WinRight] = maskLinescanVisual(Openfile)
    % TODO: add option for automatic detection
    f = figure;

    % USER CHOOSES RELEVANT AREA FOR ANALYSIS
    % show 1 frame at a time

    framenumber = 1;

    % TODO: deal with 2D/3D?
    % TODO: this is a very slow way to do this
    maxframes = length(imfinfo(Openfile));
    Iframe = imread(Openfile, framenumber);
    Iframe = imadjust(Iframe, stretchlim(Iframe, 0));
    I = imshow(Iframe);

    % Title with frame number
    title({Openfile;['frame: ', num2str(framenumber)]}, 'Interpreter', 'none');
    drawnow;

    f.KeyPressFcn = @keyPressed;

    % TODO: use drawrectangle()
    r = drawrectangle('Color','r');
    uiwait;

    WinLeft = round(r.Position(1));
    WinRight = round(r.Position(1) + r.Position(3));

    close (f);

% TODO: just set the CData for new frames;
function keyPressed(~, evt)
    if isempty(evt.Modifier)
        switch evt.Key
            case 'space' % space for keep this box
                uiresume;
            case 'f' % forward 1 frame
                framenumber = min(framenumber + 1, maxframes);
                updateFrame();
            case 'b' % back 1 frame
                framenumber = max(framenumber - 1, 1);
                updateFrame()
            case 's' % skip 10 frames forward
                framenumber = min(framenumber + 10, maxframes);
                updateFrame()
        end
    end
end

    function updateFrame()
        Iframe = imread(Openfile, framenumber);
        Iframe = imadjust(Iframe, stretchlim(Iframe, 0));
        I.CData = Iframe;
        title({Openfile;['frame: ', num2str(framenumber)]}, 'Interpreter', 'none');
        drawnow;
    end
end
