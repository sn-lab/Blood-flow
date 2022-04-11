function f = plotLinescan(I, varargin)
    % Settings
    gridColor = [1, 0, 0];
    gridAlpha = 0.5;
    gridLineWidth = 2;
    
    
    f = figure();
    
    % Plot linescan velocity trace above
    hPlotAxes = axes('Position',[0.05 0.71 0.9 0.2]);
    plot(hPlotAxes, varargin{:})
    ylabel('Velocity (mm/s)');
    hPlotAxes.XTickLabels = [];
    hPlotAxes.XGrid = 'on';
    hPlotAxes.GridColor = gridColor;
    hPlotAxes.GridAlpha = gridAlpha;
    drawnow;
%     hPlotAxes.XGridHandle.LineWidth = gridLineWidth;
    
    % Plot image (kymograph) below
    hImageAxes = axes('Position',[0.05 0.1 0.9 0.6]);
    % TODO: plot x axis in time coordinates and add label
    hImg = imagesc(hImageAxes, I');
    hImg.Interpolation = 'bilinear';
    colormap(hImageAxes, 'gray')
    hImageAxes.XGrid = 'on';
    hImageAxes.GridColor = gridColor;
    hImageAxes.GridAlpha = gridAlpha;
    drawnow;
%     hImageAxes.XGridHandle.LineWidth = gridLineWidth;

    % Use linkaxes/linkprop to link x axes
    linkaxes([hPlotAxes, hImageAxes], 'x');
end