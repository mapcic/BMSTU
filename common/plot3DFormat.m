function [  ] = plot3DFormat( Axes, shadingVal, colormapVal, viewVal)
    axes(Axes);
    Axes.NextPlot = 'add';
    
    if ( ~isempty(shadingVal) )
        shading(Axes, shadingVal);
    end
    
    if ( ~isempty(colormapVal) )
        colormap(Axes, colormapVal);
    end
 
    if ( ~isempty(viewVal) )
        view(Axes, viewVal);
    end
end