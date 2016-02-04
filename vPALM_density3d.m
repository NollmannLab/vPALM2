
%% density3d
function density=vPALM_density3d(h,x,y,z)
status(h,'r','working...');
tic
% det_thresh=str2num(get(h.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));
det_thresh=str2num(get(h.h_sliders.detect_thresh,'String'))./str2num(get(h.pixelsize,'String'));

parfor iparticle=1:size(x,1)
    n=find( sqrt( (x(iparticle)-x).^2 + (y(iparticle)-y).^2 + (z(iparticle)-z).^2) < det_thresh) ;
    density(iparticle)=size(n,1);
    
end

disp('density calculated');
status(h,'g',strcat('ok, i am done in ',num2str(toc),' s'));

