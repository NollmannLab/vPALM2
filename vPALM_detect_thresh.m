
    
function vPALM_detect_thresh(hObject, eventdata, h)
    set(h.detect_thresh_slider,'Value',str2num(get(h.detect_thresh,'String')));
    vPALM_plot_localizations_2d(hObject,eventdata,h)

