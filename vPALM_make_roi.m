%% make_roi
function vPALM_make_roi(hObject, eventdata, h)


axes(h.axes1)
h1 = imrect(h.axes1);
rectzoom = wait(h1); % wait until you double-click on the image and return the position of the rectangle

h.rectzoom=rectzoom;

set(h.h_sliders.Xmax,'String',num2str(h.rectzoom(1)+rectzoom(3)));
set(h.h_sliders.Xmin,'String',num2str(h.rectzoom(1)));
set(h.h_sliders.Ymax,'String',num2str(h.rectzoom(2)+rectzoom(4)));
set(h.h_sliders.Ymin,'String',num2str(h.rectzoom(2)));

if h.rectzoom(1)+rectzoom(3)>get(h.h_sliders.Xmax_slider,'Max')
    set(h.h_sliders.Xmax_slider,'Max',h.rectzoom(1)+rectzoom(3)+1)
end
set(h.h_sliders.Xmax_slider,'Value',h.rectzoom(1)+rectzoom(3));

if h.rectzoom(1)<get(h.h_sliders.Xmin_slider,'Min')
    set(h.h_sliders.Xmin_slider,'Min',h.rectzoom(1)-1)
end
set(h.h_sliders.Xmin_slider,'Value',h.rectzoom(1));

if h.rectzoom(2)+rectzoom(4)>get(h.h_sliders.Ymax_slider,'Max')
    set(h.h_sliders.Ymax_slider,'Max',h.rectzoom(2)+rectzoom(4)+1)
end
set(h.h_sliders.Ymax_slider,'Value',h.rectzoom(2)+rectzoom(4));

if h.rectzoom(2)<get(h.h_sliders.Ymin_slider,'Min')
    set(h.h_sliders.Ymin_slider,'Min',h.rectzoom(2)-1)
end
set(h.h_sliders.Ymin_slider,'Value',h.rectzoom(2));

setcallbacks(h);


vPALM_plot_localizations_2d(hObject,eventdata,h)
