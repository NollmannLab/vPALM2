%% 
function m2=filter_sliders(h,m)
    p=read_parameters(h);
    
    if h.parameters.fit_3dxcorr==0 & h.parameters.fit_3dgaussian==0
        
        if get(h.chb_drift_correction,'Value')==0

            r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
                    & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax);
            m2=m(r,:);   
        else %uses drift correction
            m=h.m_drift;
            r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
                    & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax);
            m2=m(r,:);   

        end
%         r1=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax);
%         r2=find(m(r1,4)>p.Imin & m(r1,4)<p.Imax);
%         r3=find(m(r2,2)>p.Xmin & m(r2,2)<p.Xmax);
%         r4=find(m(r3,3)>p.Ymin & m(r3,3)<p.Ymax);
%         m2=m(r4,:);


    elseif h.parameters.fit_3dxcorr==1 %      plots 3d xcorr localizations

        r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
            & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax ...
            & m(:,13)>p.Zmin & m(:,13)<p.Zmax);
        m2=m(r,:);

    elseif h.parameters.fit_3dgaussian==1 %      plots 3d gaussian localizations    
   
        if h.parameters.v_gauss_cal==0 % no calibation file, wil show ratio wy/wx;
            z=m(:,13)./m(:,14);
            m(:,15)=z;
            disp('plotting wx/wy as no cal is available')
        else % uses cal
%             z=calc_z_from_cal(h,m,h.parameters);
%             m(:,15)=z;
            if get(h.chb_3dcalapplied,'Value')==1
                m=h.m_3d;
                disp('using calibration')
            else
                m=h.m;
                z=m(:,13)./m(:,14);
                m(:,15)=z;
                disp('no calibration') 
            end
        end

        r=find(m(:,1)>p.Frmin & m(:,1)<p.Frmax & m(:,2)>p.Xmin & m(:,2)<p.Xmax ...
            & m(:,3)>p.Ymin & m(:,3)<p.Ymax & m(:,4)>p.Imin & m(:,4)<p.Imax ...
            & m(:,15)>p.Zmin & m(:,15)<p.Zmax);
        m2=m(r,:);

end

