% plotFormat(Axes{i}, 'Title', 'xTitle', 'yTitle', {}, [], [])

function [  ] = plotFormat( Axes, title, xTitle, yTitle, legends, xLim, yLim )
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
    
    Axes.Title.String = title;
    Axes.XLabel.String = xTitle;
    Axes.YLabel.String = yTitle;
    
    if ( ~isempty(legends) )
        legend(legends, 'Interpreter', 'latex');
    end
end