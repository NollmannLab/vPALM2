
%%    vPALM
function z=calc_z_from_cal(h,m,parameters);


        ax=h.zcal_xpars;
        ay=h.zcal_ypars;
        zgrid=linspace(parameters.cal_from,parameters.cal_to,1000);
        wx_cal=polyval(ax,zgrid);
        wy_cal=polyval(ay,zgrid);
        wx=m(:,13);
        wy=m(:,14);
      
        
        minsq2=[];
        for ievent=1:size(m,1)
            [min2(ievent),minsq2(ievent)]=min( (wx(ievent)-wx_cal).^2+(wy(ievent)-wy_cal).^2 );
        end
        z=(zgrid(minsq2)-parameters.cal_from)*str2num(get(h.stepsize,'String'));
        z=z-min(z);
%         c=z;
disp('done')

