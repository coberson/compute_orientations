function showVectorfield3D( volume,vectorfield_1,vectorfield_2 )
    %volume has to be logical matrix
    %vectorfield has to be a 4-dimensional matrix
    
    if nargin<3
        vectorfield_2=zeros(size(vectorfield_1));
    end
    figure
    [X,Y] = meshgrid(1:size(volume,2),1:size(volume,1));

    for i=1:1:size(volume,3)
        fprintf('Z: %d\n',i);
        clf
        v = double(volume>0);
        p = patch( isosurface(v,0) ); %# create isosurface patch
        isonormals(v, p);
        axis([0 size(volume,2) 0 size(volume,1) 0 size(volume,3)]);set(p, 'FaceColor','r', 'EdgeColor','none');grid on;
        alpha(0.5);lighting gouraud;
        daspect([1 1 1])
        hold on
        quiver3(X,Y,zeros(size(X))+i,vectorfield_1(:,:,i,2),vectorfield_1(:,:,i,1),vectorfield_1(:,:,i,3),1);rotate3d on;
        quiver3(X,Y,zeros(size(X))+i,vectorfield_2(:,:,i,2),vectorfield_2(:,:,i,1),vectorfield_2(:,:,i,3),1);rotate3d on;
        pause %caution: I permuted the two components
    end
end

