

%% selectbeads
function vPALM_selectbeads(hObject, eventdata, h)

% h=vPALM2_Beads_detection(h);
[h,PointsToRemove]=vPALM2_Beads_detection_JB_v1(h);

AnalysisMethods = get(h.load_options, 'Value'); % Return the method used to detect the single particle fluorescent events
Applied3DCalibration = get(h.chb_cal_3d, 'Value');

if AnalysisMethods && Applied3DCalibration
    m = h.m_3d;
    m(:,2) = m(:,2) - h.Ref_position(m(:,1),1);
    m(:,3) = m(:,3) - h.Ref_position(m(:,1),2);
    m(:,15) = m(:,15) - h.Ref_position(m(:,1),3);
else
    m = h.m;
    m(:,2) = m(:,2) - h.Ref_position(m(:,1),1);
    m(:,3) = m(:,3) - h.Ref_position(m(:,1),2);
end

if get(h.chb_Remove_beads, 'Value')
    m(PointsToRemove,:) = [];
end

if Applied3DCalibration
    h.m_3d_drift = m;
else
    h.m_drift = m;
end


% Update handles structure
setcallbacks(h)


% m2=h.m;
% m2(:,2)=m2(:,2)-h.Ref_position(m2(:,1),1);
% m2(:,3)=m2(:,3)-h.Ref_position(m2(:,1),2);
% h.m_drift=m2;
% 
% set(h.chb_drift_correction,'Value',1);
% setcallbacks(h);
% disp('done');
% 
% % make_roi(hObject, eventdata, h);
% % m=filter_sliders(h,h.m);
% % 
% % axes(h.axes2);
% % 
% % I=m(:,4);
% % maxI=max(I);
% % 
% % select=find(I>0);%get(h.showbeads_slider,'Value')*maxI);
% % t=m(select,1);
% % x=m(select,2);
% % y=m(select,3);
% % meanx=mean(x);
% % meany=mean(y);
% % 
% % xm=meanfilter(100,x)'+meanx;
% % ym=meanfilter(100,y)'+meany;
% % 
% % 
% % 
% % %makes interpolation
% % frmin=min(h.m(:,1));
% % frmax=max(h.m(:,1));
% % frrange=[frmin:1:frmax+1];
% % h.tm=frrange;
% % xm=interp1(t,xm,frrange);
% % ym=interp1(t,ym,frrange);
% % h.xm(frrange)=xm';
% % h.ym(frrange)=ym';
% % 
% % x2=x-h.xm(t)';
% % y2=y-h.ym(t)';
% % 
% % % for i=1:size(x,1)
% % %     x2(:)=x(:)-h.xm(find(h.tm==t(:)));
% % %     y2(:)=y(:)-h.ym(find(h.tm==t(:)));
% % % end
% % 
% % hold off,
% % plot(x2,y2,'.b')
% % 
% % x2=x-h.xm(t)';
% % y2=y-h.ym(t)';