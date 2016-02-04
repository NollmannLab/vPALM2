


function vPALM_apply_3dcal(hObject, eventdata, h)

if get(h.chb_cal_3d,'Value')==1
    h.m_3d=vPALM2_calc_3d(h);

    set(h.chb_3dcalapplied,'Value',1);
else
    disp('No calibration was done!');
end

% Update handles structure
setcallbacks(h)