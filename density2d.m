%% density2d
function density=density2d(h,x,y)
status(h,'r','working...')
% get(h.h_sliders.Xmax_slider,'Max')
% det_thresh=str2num(get(h.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));
det_thresh=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

parfor iparticle=1:size(x,1)
    n= find( sqrt( (x(iparticle)-x).^2 + (y(iparticle)-y).^2 ) < det_thresh);
    density(iparticle)=size(n,1 );
    
end

disp('density calculated');
    status(h,'g','ok, i am done')
