function m=filter_3dgauss(m0,parameters,h)

if get(h.load_options,'Value')==6

    ellip=m0(:,11);
    
    mean_w=m0(:,13)/h.pixelsize;

else
    
    ellip=m0(:,13)./m0(:,14);
    
    mean_w=sqrt(m0(:,13).*m0(:,14));
    
end

    keep_values=find(ellip>parameters.min_ellip & ellip<parameters.max_ellip & mean_w>parameters.min_w & mean_w< parameters.max_w);

    m=m0(keep_values,:);