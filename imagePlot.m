function imagePlot(I, varargin)
    f = figure;
    
    % 
    subplot(2,1,1);
    plot(varargin)
    
    %
    subplot(2,1,2);
    imshow(I, []);
    
    % Use linkaxes/linkprop to link x axes
end