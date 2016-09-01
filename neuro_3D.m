%% test 3D orientation on the neuroblasts of sample seg_Tm2
clear all
close all
clc

FileNameWP='/share/numeric/imaging_data/neuroblasts/seg_mat_stacks_2011/seg_Tm2';
FileName=FileNameWP;
path=strfind(FileName,'/');

if length(path)>0
    FileName=FileName(path(length(path))+1:end);
end

fileopen=strcat(FileNameWP,'.mat');
m=load(fileopen);
gm=m.gm; %matrix containing the green intensities
rm=m.rm; %matrix containing the red intensities
L=m.newdmatL; %CHANGE output of program that we get all these three matrices in newdmatL!
cpnts=unique(L);

for i=1:max(cpnts)
    if i~=2
        fprintf('Component number: %u\n',i);
        comp=(L==i).*gm;
        
        radius=10;
        
        % define minimal box for the component
        [nx,ny,nz]=size(comp);
        [r,c]=find(comp);
        [c1,c2]=ind2sub([ny,nz],c);
        minx=min(r);
        maxx=max(r);
        miny=min(c1);
        maxy=max(c1);
        minz=min(c2);
        maxz=max(c2);
        comp=comp(minx:maxx,miny:maxy,minz:maxz);
        [nx,ny,nz]=size(comp);
        
        fprintf('Dimensions: %u %u %u \n',nx,ny,nz);
        [comp,v1]=compute_orientations_3D(comp,radius);
        sizev=size(v1);
        fprintf('Dimensions: %u %u %u\n',sizev(1),sizev(2),sizev(3));
    end
end
