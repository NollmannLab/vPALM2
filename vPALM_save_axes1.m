%% save_axes1
function vPALM_save_axes1(hObject, eventdata, h)

axes(h.axes1)
% export_fig test3.png
F=getframe(h.axes1);               %select axes in GUI
figure();                                          %new figure
image(F.cdata);  %show selected axes in new figure
saveas(gcf, strcat(h.fullFileName,'_axes1.fig'), 'fig');                    %save figure
close(gcf); %and close it

% 
% axes(h.axes1)
% h1=gcf;`
% h2=figure;
% objects=allchild(h1);
% copyobj(get(h1,'children'),h2);

%% save_axes2
% function save_axes2(hObject, eventdata, h)
% 
% axes(h.axes2)
% % export_fig test3.png
% F=getframe(h.axes2);               %select axes in GUI
% figure();                                          %new figure
% image(F.cdata);  %show selected axes in new figure
% saveas(gcf, strcat(h.fullFileName,'_axes2.png'), 'png');                    %save figure
% close(gcf);                                       %and close it
