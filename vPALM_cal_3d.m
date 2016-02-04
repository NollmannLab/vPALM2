

%%
function vPALM_cal_3d(hObject, eventdata, h)

tic
[h.wxPoly,h.wyPoly,h.ccPolyWtoZ]=vPALM2_cal3D(h);

set(h.chb_cal_3d,'Value',1);
    
setcallbacks(h);
status(h,'g',strcat('3d data plotted in :',num2str(toc,2),' s'));
toc
disp('done') ;