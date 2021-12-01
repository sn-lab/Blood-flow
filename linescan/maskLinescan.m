function [WinLeft, WinRight] = maskLinescan(I)
% TODO: add option for automatic detection
f = figure;

% USER CHOOSES RELEVANT AREA FOR ANALYSIS
% show 1 frame at a time

framenumber = 1;

% TODO: deal with 2D/3D?
imshow(I, []);
title({fname;['frame:', num2str(framenumber)]});

% get coordinates of rrbox
xlabel('Select region of interest');

f.KeyPressFcn = @keyPressed;

% TODO: use drawrectangle()
r = drawrectangle('Color','r');
uiwait;

WinLeft = r.Position(1);
WinRight = r.Position(1) + r.Position(3);

% DO this is 
close (f);

end

% TODO: just set the CData for new frames;
function keyPressed(src, evt)
    if isempty(evt.Modifier)
        switch evt.Key
            case ' ' % space for keep this box
                uiresume;
            case 'f' % forward 1 frame
                if framenumber < maxframes
                    framenumber = framenumber +1;
                    showlines = imread(Openfile, framenumber);
                else
                    beep;
                    disp('no more frames')
                end
                
                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber)]});

            case 'b' % back 1 frame
                if framenumber > 1
                    framenumber = framenumber - 1;
                    showlines = imread(Openfile, framenumber);
                else
                    beep;

                end

                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber)]});

            case 's' % skip 10 frames forward
                if framenumber + 10 <= maxframes
                    framenumber = framenumber + 10;
                    showlines = imread(Openfile, framenumber);
                else
                    beep;
                end

                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber)]});

            case 'n'
                clf;
                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber)]});
        end
    end
end

