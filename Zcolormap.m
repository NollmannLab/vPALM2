%%
function Zcolormap(hObject, eventdata, h)

v=get(h.Zcolormap,'Value');

switch v
    case 1
        colormap jet
    case 2
        colormap hot
    case 3
        colormap HSV
    case 4
        colormap gray
end

% plot_localizations(hObject,eventdata,h);
