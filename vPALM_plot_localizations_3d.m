%% plot_localizations 3d
function h=vPALM_plot_localizations_3d(hObject,eventdata,h)
    status(h,'r','working...')

if get(h.chb_data_loaded,'Value')==1
    
    m0=h.m;
    parameters=h.parameters;
    v=get(h.Zcolormap,'Value');
    switch v
        case 1
            parameters.colormap='jet';
        case 2
            parameters.colormap='hot';
        case 3
            parameters.colormap='hsv';
        case 4
            parameters.colormap='gray';
        case 5
            parameters.colormap='cool';            
    end
    
    if parameters.fit_2dgaussian==1 | (parameters.fit_3dgaussian==0 & parameters.fit_3dxcorr==0)
        vPALM_plot_m_2d(h,m0,parameters)
    elseif (parameters.fit_3dgaussian==1 | parameters.fit_3dxcorr==1)
        vPALM_plot_m_3d(h,m0,parameters);
    end

else
    status(h,'r','No file loaded!')
end

% Update handles structure
setcallbacks(h)

