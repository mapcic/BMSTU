function [  ] = plotFormat( Axes, title, xTitle, yTitle, legends, xLim, yLim, textSize )
    axes(Axes);
    grid on;
    Axes.NextPlot = 'add';
    
    Axes.Title.Interpreter = 'latex';
    Axes.XLabel.Interpreter = 'latex';
    Axes.YLabel.Interpreter = 'latex';
    
    if ( ~isempty(xLim) )
        Axes.XLim = xLim;
    end
    
    if ( ~isempty(yLim) )
        Axes.YLim = yLim;
    end
    
    if ( ~isempty(textSize) )
        Axes.XLabel.FontSize = textSize;
        Axes.YLabel.FontSize = textSize;
        Axes.Title.FontSize = textSize + 4;
    end

    Axes.Title.String = title;
    Axes.XLabel.String = xTitle;
    Axes.YLabel.String = yTitle;
    
    if ( ~isempty(legends) )
        L = legend(legends, 'Interpreter', 'latex');
        if ( ~isempty(textSize) )
            set(L, 'FontSize', textSize);
        end
    end
end