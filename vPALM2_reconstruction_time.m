function vPALM2_reconstruction_time(x,y,t,s,ps,px)

% used to put a control particle and make sure it shows up in all channels
% aligned. removed for normal operation

% meanX=mean(x);
% meanY=mean(y);
% x(1:30)= ones(1,30)*meanX;
% y(1:30)= ones(1,30)*meanY;
% x(701:730)= ones(1,30)*meanX;
% y(701:730)= ones(1,30)*meanY;
% x(1601:1630)= ones(1,30)*meanX;
% y(1601:1630)= ones(1,30)*meanY;

% puts localizations in a matrix and calculates zooming factor
loc(:,1) = x-min(x);
loc(:,2) = y-min(y);
superzoom = ceil(px/ps);
szx = superzoom * ceil(max(loc(:,1)));
szy = superzoom * ceil(max(loc(:,2)));
adjust_threshold=.99;
%% calculates 2d histogram based on localizations
figure('Position',[600 200 600 500],'Color',[1 1 1],'name','Full reconstruction');
[hh,centers]=hist3(loc,[szx,szy]);
h = fspecial('gaussian', [5*ceil(s/ps) 5*ceil(s/ps)], ceil(s/ps));
PALM_image_full=(imfilter(hh,h));

    
    PALM_image_full=imadjust( PALM_image_full,[0 adjust_threshold],[0 1]);
    ax=imagesc(flipud(rot90(PALM_image_full)));
    set(gca,'YDir','normal')
    colormap hot
    axis equal    
    colormap('hot');
    h1 = colorbar;
    ylabel(h1, 'Density','fontsize',12)

    xlabel('X, px');
    ylabel('Y, px');

%% calculates series of images (9 as a default) and makes a reconstruction for each stack
    Nfigs=9;
    n1=ceil(sqrt(Nfigs));
    figure('Position',[800 200 1200 1000],'Color',[1 1 1],'name','Partial reconstructions by time');

    Nlist=[1:ceil(max(t)/Nfigs):max(t)];
    PALM_image=[];
    for i=1:length(Nlist)-1
        subplot(n1,n1,i)
        loc_t=[];
        range=find(t>ceil(Nlist(i)) & t<ceil(Nlist(i+1)));

%         loc_t(:,1) = x(Nlist(i):Nlist(i+1))-min(x);
%         loc_t(:,2) = y(Nlist(i):Nlist(i+1))-min(y);
        loc_t(:,1) = x(range)-min(x);
        loc_t(:,2) = y(range)-min(y);

        hh=hist3(loc_t,centers);
        PALM_image(:,:,i)=flipud(rot90((imfilter(hh,h))));
        PALM_image(:,:,i)=imadjust( PALM_image(:,:,i),[0 adjust_threshold],[0 1]);
        ax=imagesc(PALM_image(:,:,i));
        set(gca,'YDir','normal')
%         axis equal    
        i
        xlabel('X, px');
        ylabel('Y, px');
        title(strcat('Frame=',num2str(ceil(Nlist(i))),'-',num2str(ceil(Nlist(i+1)))));
    end
    subplot(n1,n1,i+1)
    imagesc(flipud(rot90(PALM_image_full)));
    set(gca,'YDir','normal')
    colormap('hot');
    xlabel('X, px');
    ylabel('Y, px');

%% splits all localizations in 3 and makes an RGB reconstruction to see evolution of structure over time
    figure('Position',[600 200 600 500],'Color',[1 1 1],'name','RGB image of time evolution');

    Nfigs=3;
    Nlist=[1:ceil(max(t)/Nfigs):max(t)];
    PALM_image_RGB=[];
    time_total=max(t);
        range_first_third=find(t<ceil(time_total/3));
        range_second_third=find(t>ceil(time_total/3) & t<ceil(2*time_total/3));
        range_last_third=find(t>ceil(2*time_total/3));

        loc_t1(:,1) = x(range_first_third)-min(x);
        loc_t1(:,2) = y(range_first_third)-min(y);        
        PALM_image_RGB(:,:,1)=flipud(rot90((imfilter(hist3(loc_t1,centers),h))));
%         PALM_image_RGB(:,:,1)=flipud(rot90((hist3(loc_t1,centers))));
        PALM_image_RGB(:,:,1)=imadjust(PALM_image_RGB(:,:,1),[0 adjust_threshold],[0 1]);

        loc_t2(:,1) = x(range_second_third)-min(x);
        loc_t2(:,2) = y(range_second_third)-min(y);
        %         PALM_image_RGB(:,:,2)=flipud(rot90((hist3(loc_t2,centers))));
        PALM_image_RGB(:,:,2)=flipud(rot90((imfilter(hist3(loc_t2,centers),h))));
        PALM_image_RGB(:,:,2)=imadjust(PALM_image_RGB(:,:,2),[0 adjust_threshold],[0 1]);


        loc_t3(:,1) = x(range_last_third)-min(x);
        loc_t3(:,2) = y(range_last_third)-min(y);
%         PALM_image_RGB(:,:,3)=flipud(rot90((hist3(loc_t3,centers))));
        PALM_image_RGB(:,:,3)=flipud(rot90((imfilter(hist3(loc_t3,centers),h))));
        PALM_image_RGB(:,:,3)=imadjust(PALM_image_RGB(:,:,3),[0 adjust_threshold],[0 1]);
    
    ax=imagesc((PALM_image_RGB))
    set(gca,'YDir','normal')
    xlabel('X, px');
    ylabel('Y, px');
    axis equal
    h1 = colorbar;
    ylabel(h1, 'Time','fontsize',12)

    
    
end
