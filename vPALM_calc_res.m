
%%
function vPALM_calc_res(hObject, eventdata, h)

status(h,'g','Calculating image resolution')
% disp(strcat('Data saved to: ',filename))

h.m2=filter_sliders(h,h.m);


ps=str2num(get(h.stepsize_reconstruction,'String'));
px = str2num(get(h.pixelsize,'String'));

x=h.m2(:,2);
y=h.m2(:,3);
t=h.m2(:,1);

coords(:,1) = x-min(x);
coords(:,2) = y-min(y);
coords(:,3) = t;

superzoom = ceil(px/ps);
szx = superzoom * ceil(max(coords(:,1)));
szy = superzoom * ceil(max(coords(:,2)));

im = binlocalizations(coords, szx, szy, superzoom);
hh=dipshow(im);
dipmapping(hh,[0 5],'colormap',hot)

fprintf('\n -- computing FIRE --\n')
[fire_value, ~, fireH, fireL] = postoresolution(coords, szx, superzoom); % in super-resolution pixels
fprintf('FIRE value %2.1f +- %2.2f [px]\n', fire_value, (fireL-fireH)/2);
fprintf('FIRE value %2.1f +- %2.2f [nm]\n', fire_value*px/superzoom, (fireL-fireH)/2*px/superzoom);


%%
%% compute FRC curve
fprintf('\n -- computing FRC curve--\n')
[~,frc_curve] = postoresolution(coords, szx, superzoom); 
figure;
qmax = 0.5/(px/superzoom);
plot(linspace(0,qmax*sqrt(2), length(frc_curve)), frc_curve,'-')
xlim([0,qmax])
hold on
plot([0 qmax],[1/7 1/7],'r-');
plot([0 qmax],[0 0],'k--'); hold off
xlabel('spatial frequency (nm^{-1})')
ylabel('FRC')
title('Fourier Ring Correlation curve')


status(h,'g',strcat('FIRE value :',num2str(fire_value*px/superzoom),' +- ',num2str((fireL-fireH)/2*px/superzoom))) 
