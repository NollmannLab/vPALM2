%% read_parameters
function p=read_parameters(h)
    p.Zmax=str2num(get(h.h_sliders.Zmax,'String'));
    p.Zmin=str2num(get(h.h_sliders.Zmin,'String'));
    p.Xmax=str2num(get(h.h_sliders.Xmax,'String'));
    p.Xmin=str2num(get(h.h_sliders.Xmin,'String'));
    p.Ymax=str2num(get(h.h_sliders.Ymax,'String'));
    p.Ymin=str2num(get(h.h_sliders.Ymin,'String'));
    p.Frmax=str2num(get(h.h_sliders.Frmax,'String'));
    p.Frmin=str2num(get(h.h_sliders.Frmin,'String'));
    p.Imax=str2num(get(h.h_sliders.Imax,'String'));
    p.Imin=str2num(get(h.h_sliders.Imin,'String'));

%     set(h.h_sliders.detect_thresh,'String','-10');

