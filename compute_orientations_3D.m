function [f,v1,v2] = compute_orientations_3D(f,radius,put_min_box)
%COMPUTE_ORIENTATIONS_3D This function computes the local structure
% orientation for a given density 3D image based on the structure tensor
% It is a complete version: In case of double orientation vector, both are given.
% v1 and v2 are chosen to be unitary.
% The radius is important for the inner product defining the structure
% tensor. It is the only parameter to adjust and it defines how big a neighbourhood we consider for every voxel.
% If put_min_box is 1, the image is cropped before computation. 
% For the computation, the image is included in a
% bigger one in order to compute the inner product with the given radius even for the boundary voxels.
%
%   Author: Chantal Oberson Ausoni, Workgroup Imaging, University of Muenster.
%   30.3.2012
%
% Needs function 'fspecial3' (Damien Garcia, in '~/matlab').
% Needs function 'showVectorfield3D(f,v1,v2)' (Hendrik
% Dirks, in '~/matlab') for 3D representation of the horizontal layers.
% Needs my function 'ellipsoid_gf' for one of the examples.
% Reference for eigenvalue computation: Oliver K. Smith Eigenvalues of a symmetric 3x3  matrix
% As f can be modified (put in a minimal box), it can be get as an output
% parameter. 
% A call 
% >>[f,v1,v2] = compute_orientations_3D; 
% executes the function for the simple example of a
% brick shape, an ellipsoid or two orthogonal ellipsoids.
% To get a 3D representation of the vectorfields: 
% >>[f,v1,v2] = s_compute_orientations_3D(f,radius);
% >>showVectorfield3D(f,v1,v2);

%% simple example
if nargin<3
    put_min_box=0;
end

% if nargin == 0 %brick example
%     f=zeros(100,120,98); 
%     f(45:54,20:100,24:54)=1; 
% end

% if nargin == 0 %ellipsoids' cross example
%     f = zeros(101,121,99); 
%     [x,y,z] = ndgrid(-50:50,-60:60,-49:49);
%     I = (x.*x/5^2+y.*y/30^2+z.*z/5)<=1; %first ellipsoid
%     f(I) = 1;
%     J = (x.*x/30^2+y.*y/5^2+z.*z/5)<=1; %second ellipsoid (orthogonal to the first)
%     f(J) = 1;
%     clear I J x y z
%     radius = 10;
%     put_min_box = 1;
% end

if nargin==0 %ellipsoid threads' example (sausages' pack)
    f= zeros(101,121,99); 
    for i=-2:2
        I=ellipsoid_gf(-i*5,i*5,i*10,5,40,5,3*pi/4,pi/12,50,60,49);
        f=f+I;
    end
    I=ellipsoid_gf(5,5,0,5,40,5,3*pi/4,7*pi/12,50,60,49);
    f=f+I;
    f=sign(f);
    radius = 10;
    put_min_box = 1;
    clear I
end


figure
v = double(f>0);
p = patch( isosurface(v,0) );
isonormals(v, p);
set(p, 'FaceColor','r', 'EdgeColor','none');
axis([0 size(f,2) 0 size(f,1) 0 size(f,3)]);
grid on;
alpha(0.5);lighting gouraud;
daspect([1 1 1])
title('Object');
clear p v

tic;


%% define minimal box
if put_min_box==1
    [~,ny,nz]=size(f);
    [r,c]=find(f);
    [c1,c2]=ind2sub([ny,nz],c);
    minx=min(r);
    maxx=max(r);
    miny=min(c1);
    maxy=max(c1);
    minz=min(c2);
    maxz=max(c2);
    f=f(minx:maxx,miny:maxy,minz:maxz);
    
    clear r c c1 c2 minx miny maxx maxy minz maxz
end 


%% include image in a bigger box to allow computation of the points on the boundaries of the box
[nx,ny,nz]=size(f);
f=cat(3,zeros(2*radius+nx,2*radius+ny,radius),cat(2,zeros(2*radius+nx,radius,nz),cat(1,zeros(radius,ny,nz),f,zeros(radius,ny,nz)),zeros(2*radius+nx,radius,nz)),zeros(2*radius+nx,2*radius+ny,radius));


%% compute derivatives
[nx ny nz] = size(f);
f=double(f);
filtx=zeros(3,3,3);
filtx(:,2,2)=[-1;0;1];
fx=imfilter(f,filtx,'symmetric','same');
filty=zeros(3,3,3);
filty(2,:,2)=[-1;0;1];
fy=imfilter(f,filty,'symmetric','same');
filtz=zeros(3,3,3);
filtz(2,2,:)=[-1;0;1];
fz=imfilter(f,filtz,'symmetric','same'); 
%fz=1/3.7939*fz; %scale with the size of the voxels????????????????????
clear filtx filty filtz


%% define weighted inner product (on a non spherical ellipsoid, because of the size of the voxels: 0.263, 0.263, 1)
% seems to give a better invariance through rotation than 'average'
radius_z=radius; 
%radius_z=ceil(radius/3.7939);
kernel = fspecial3('ellipsoid', [2*radius+1,2*radius+1,2*radius_z+1]);  %caution: is constant on the ball (no interpolation for the boundary!!!)
w_inner_prod = @(f,g) convn(f.*g,kernel,'same');


%% define structure tensor
a = w_inner_prod(fx,fx); 
b = w_inner_prod(fy,fy);
c = w_inner_prod(fz,fz);
ab = w_inner_prod(fx,fy);
ac = w_inner_prod(fx,fz);
bc = w_inner_prod(fy,fz);


%% define useful parameters
dst=2.*ab.*ac.*bc+a.*b.*c-a.*bc.^2-b.*ac.^2-c.*ab.^2; %determinant of the structure tensor
m=1/3.*(a+b+c); 
q = (m.*(ab.^2+ac.^2+bc.^2-a.*c-b.*c-a.*b)+m.^2.*(a+b+c)-m.^3+dst) / 2;
%q=q/2;
p = (2*ab.^2 +2*ac.^2 +2*bc.^2+(a-m).^2 +(b-m).^2 +(c-m).^2) / 6;
%p=1/6*p;
trexp=2*ab.^2+2*ac.^2+2*bc.^2-2*a.*b-2*c.*a-2*b.*c; %for the multiplicity, further


%% compute eigenvalues
newp = (p~=0).*(p.^1.5)+(p==0).*ones(nx,ny,nz);
phi  = (p~=0).*(abs(q)<=abs(p.^1.5)).*(1/3).*acos(q./newp)+(p~=0).*(abs(q)>abs(p.^1.5)).*(q.*p<0).*pi/3.*ones(size(a)); 
%the case p=0 implies the matrix is m*identity, if the quotient is not
%between -1 and 1, the arccos cannot be taken. If the quotient is more than 1
%use 0 instead, if it is less than -1 use pi/3 instead! 

lambda = m-sqrt(p).*cos(phi)-sqrt(p).*sqrt(3).*sin(phi);
lambda=(dst~=0).*lambda; %dst>epsilon?

clear newp phi 


%% write cross products to obtain one or two vectors perpendicular to the image of mat-lambda*id
epsilon=10^(-15); %parameter for the precision

norm1=(a-lambda).^2+ab.^2+ac.^2;
norm2=ab.^2+(b-lambda).^2+bc.^2;
norm3=ac.^2+bc.^2+(c-lambda).^2;

% columns of A-lambda*I, sort them to have the zero ones at the end
colspace1(:,:,:,1) = (norm1>epsilon).*(a-lambda)+(norm1<epsilon).*(norm2>epsilon).*ab+...
    (norm1<epsilon).*(norm2<epsilon).*(norm3>epsilon).*ac;
colspace1(:,:,:,2) = (norm1>epsilon).*ab+(norm1<epsilon).*(norm2>epsilon).*(b-lambda)+...
    (norm1<epsilon).*(norm2<epsilon).*(norm3>epsilon).*bc;
colspace1(:,:,:,3) = (norm1>epsilon).*ac+(norm1<epsilon).*(norm2>epsilon).*bc+...
    (norm1<epsilon).*(norm2<epsilon).*(norm3>epsilon).*(c-lambda);
ncs1=colspace1(:,:,:,1).^2+colspace1(:,:,:,2).^2+colspace1(:,:,:,3).^2;%squared norm!

ncs1=(ncs1>0).*ncs1+(ncs1==0).*ones(nx,ny,nz);

coeff=1./sqrt(ncs1);
colspace1=cat(4,coeff,coeff,coeff).*colspace1;

clear coeff
ncs1=colspace1(:,:,:,1).^2+colspace1(:,:,:,2).^2+colspace1(:,:,:,3).^2;%squared norm 0 or 1

colspace2(:,:,:,1) = (norm1>epsilon).*(norm2>epsilon).*ab+...
    (norm1<epsilon).*(norm2>epsilon).*ac+...
    (norm1>epsilon).*(norm2<epsilon).*(norm3>epsilon).*ac;
colspace2(:,:,:,2) = (norm1>epsilon).*(norm2>epsilon).*(b-lambda)+...
    (norm1<epsilon).*(norm2>epsilon).*bc+...
    (norm1>epsilon).*(norm2<epsilon).*(norm3>epsilon).*bc;
colspace2(:,:,:,3) = (norm1>epsilon).*(norm2>epsilon).*bc+...
    (norm1<epsilon).*(norm2>epsilon).*(c-lambda)+...
    (norm1>epsilon).*(norm2<epsilon).*(norm3>epsilon).*(c-lambda);
ncs2=colspace2(:,:,:,1).^2+colspace2(:,:,:,2).^2+colspace2(:,:,:,3).^2;%squared norm!

colspace3(:,:,:,1) = (norm1>epsilon).*(norm2>epsilon).*(norm3>epsilon).*ac;
colspace3(:,:,:,2) = (norm1>epsilon).*(norm2>epsilon).*(norm3>epsilon).*bc;
colspace3(:,:,:,3) = (norm1>epsilon).*(norm2>epsilon).*(norm3>epsilon).*(c-lambda);
ncs3=colspace3(:,:,:,1).^2+colspace3(:,:,:,2).^2+colspace3(:,:,:,3).^2;%squared norm!

clear norm1 norm2 norm3
clear a b c ab ac bc lambda 

%cross products in case lambda has multiplicity one
normprod=(ncs1>epsilon).*(ncs2>epsilon).*ncs1.*ncs2+(ncs1<epsilon).*ones(nx,ny,nz)+(ncs2<epsilon).*(ncs1>epsilon).*ones(nx,ny,nz);
normprod=1./sqrt(normprod);

cross1(:,:,:,1) = normprod.*(colspace1(:,:,:,2).*colspace2(:,:,:,3)-colspace1(:,:,:,3).*colspace2(:,:,:,2));
cross1(:,:,:,2) = normprod.*(colspace1(:,:,:,3).*colspace2(:,:,:,1)-colspace1(:,:,:,1).*colspace2(:,:,:,3));
cross1(:,:,:,3) = normprod.*(colspace1(:,:,:,1).*colspace2(:,:,:,2)-colspace1(:,:,:,2).*colspace2(:,:,:,1));
normc1=cross1(:,:,:,1).^2+cross1(:,:,:,2).^2+cross1(:,:,:,3).^2;%squared norm!

normprod=(ncs1>epsilon).*(ncs3>epsilon).*ncs1.*ncs3+(ncs1<epsilon).*ones(nx,ny,nz)+(ncs3<epsilon).*(ncs1>epsilon).*ones(nx,ny,nz);
normprod=1./sqrt(normprod);
clear colspace2 ncs2

cross2(:,:,:,1) = normprod.*(colspace1(:,:,:,2).*colspace3(:,:,:,3)-colspace1(:,:,:,3).*colspace3(:,:,:,2));
cross2(:,:,:,2) = normprod.*(colspace1(:,:,:,3).*colspace3(:,:,:,1)-colspace1(:,:,:,1).*colspace3(:,:,:,3));
cross2(:,:,:,3) = normprod.*(colspace1(:,:,:,1).*colspace3(:,:,:,2)-colspace1(:,:,:,2).*colspace3(:,:,:,1));
normc2=cross2(:,:,:,1).^2+cross2(:,:,:,2).^2+cross2(:,:,:,3).^2;%squared norm!

cross1(:,:,:,1)=(normc1>epsilon).*cross1(:,:,:,1)+(normc1<epsilon).*(normc2>epsilon).*cross2(:,:,:,1);
cross1(:,:,:,2)=(normc1>epsilon).*cross1(:,:,:,2)+(normc1<epsilon).*(normc2>epsilon).*cross2(:,:,:,2);
cross1(:,:,:,3)=(normc1>epsilon).*cross1(:,:,:,3)+(normc1<epsilon).*(normc2>epsilon).*cross2(:,:,:,3);
normc1=cross1(:,:,:,1).^2+cross1(:,:,:,2).^2+cross1(:,:,:,3).^2;%squared norm!
clear cross2 colspace3 ncs3

%cross products in case lambda has multiplicity two
normprod=(ncs1>epsilon).*3.*ncs1+(ncs1<epsilon).*ones(nx,ny,nz);
normprod=1./sqrt(normprod);

cross3(:,:,:,1)=normprod.*(colspace1(:,:,:,2)-colspace1(:,:,:,3)); %cross product the column of colspace with (1,1,1)
cross3(:,:,:,2)=normprod.*(colspace1(:,:,:,3)-colspace1(:,:,:,1));
cross3(:,:,:,3)=normprod.*(colspace1(:,:,:,1)-colspace1(:,:,:,2));
normc3=cross3(:,:,:,1).^2+cross3(:,:,:,2).^2+cross3(:,:,:,3).^2;%squared norm!

normprod=(ncs1>epsilon).*ncs1+(ncs1<epsilon).*ones(nx,ny,nz);
normprod=1./sqrt(normprod);

cross4(:,:,:,1)=zeros(size(colspace1(:,:,:,1))); %in case cospace is colinear with (1,1,1), cross product of colspace with (1,0,0)
cross4(:,:,:,2)=normprod.*colspace1(:,:,:,3);
cross4(:,:,:,3)=-normprod.*colspace1(:,:,:,2);
normc4=cross4(:,:,:,1).^2+cross4(:,:,:,2).^2+cross4(:,:,:,3).^2;%squared norm!


cross3(:,:,:,1)=(normc3>epsilon).*cross3(:,:,:,1)+(normc3<epsilon).*(normc4>epsilon).*cross4(:,:,:,1);
cross3(:,:,:,2)=(normc3>epsilon).*cross3(:,:,:,2)+(normc3<epsilon).*(normc4>epsilon).*cross4(:,:,:,2);
cross3(:,:,:,3)=(normc3>epsilon).*cross3(:,:,:,3)+(normc3<epsilon).*(normc4>epsilon).*cross4(:,:,:,3);
normc3=cross3(:,:,:,1).^2+cross3(:,:,:,2).^2+cross3(:,:,:,3).^2;%squared norm!
clear cross4 normc4


%% First eigenvector for the minimal eigenvalue
v1 (:,:,:,1) = (normc1>=epsilon).*cross1(:,:,:,1)+...
    (normc1<epsilon).*(normc3>epsilon).*cross3(:,:,:,1);
v1 (:,:,:,2) = (normc1>=epsilon).*cross1(:,:,:,2)+...
    (normc1<epsilon).*(normc3>epsilon).*cross3(:,:,:,2);
v1 (:,:,:,3) = (normc1>=epsilon).*cross1(:,:,:,3)+...
    (normc1<epsilon).*(normc3>epsilon).*cross3(:,:,:,3);
n1=v1(:,:,:,1).^2+v1(:,:,:,2).^2+v1(:,:,:,3).^2;%squared norm!

n1=(n1>0).*n1+(n1==0).*ones(nx,ny,nz);
coeff=1./sqrt(n1);
v1=cat(4,coeff,coeff,coeff).*v1;
n1=v1(:,:,:,1).^2+v1(:,:,:,2).^2+v1(:,:,:,3).^2;%squared norm!
clear coeff cross1 


%% Second eigenvector for the minimal eigenvalue (in case lambda has multiplicity 2)
normprod=(ncs1>epsilon).*(normc3>epsilon).*ncs1.*normc3+(ncs1<epsilon).*ones(nx,ny,nz)+(normc3<epsilon).*(ncs1>epsilon).*ones(nx,ny,nz);
normprod=1./sqrt(normprod);
v2(:,:,:,1) = (normc1<epsilon).*(normc2<epsilon).*normprod.*(colspace1(:,:,:,2).*v1(:,:,:,3)-colspace1(:,:,:,3).*v1(:,:,:,2));
v2(:,:,:,2) = (normc1<epsilon).*(normc2<epsilon).*normprod.*(colspace1(:,:,:,3).*v1(:,:,:,1)-colspace1(:,:,:,1).*v1(:,:,:,3));
v2(:,:,:,3) = (normc1<epsilon).*(normc2<epsilon).*normprod.*(colspace1(:,:,:,1).*v1(:,:,:,2)-colspace1(:,:,:,2).*v1(:,:,:,1));
n2=v2(:,:,:,1).^2+v2(:,:,:,2).^2+v2(:,:,:,3).^2;%squared norm!

n2=(n2>0).*n2+(n2==0).*ones(nx,ny,nz);
coeff=1./sqrt(n2);
v2=cat(4,coeff,coeff,coeff).*v2;

clear coeff
n2=v2(:,:,:,1).^2+v2(:,:,:,2).^2+v2(:,:,:,3).^2;%squared norm!


clear cross3 colspace1 normprod ncs1 normc1 normc2 normc3

fprintf('Maximal norm of v1: %d\n',max(n1(:)));
fprintf('Maximal norm of v2: %d\n',max(n2(:)));


%% Cancel vectors on the boundary
%redefine weighted inner product with smaller radius
kernel = fspecial3('ellipsoid', 5);%ADJUST TO HAVE ORIENTATION FOR MORE VOXELS AROUND THE OBJECT
w_inner_prod = @(f,g) convn(f.*g,kernel,'same');
tr=f+w_inner_prod(fx,fx)+w_inner_prod(fy,fy)+w_inner_prod(fz,fz)+w_inner_prod(f,f);
coeff=(tr>0); %cancel vectors when on boundary
coeff=cat(4,coeff,coeff,coeff);
v1=coeff.*v1;
v2=coeff.*v2;
n1=v1(:,:,:,1).^2+v1(:,:,:,2).^2+v1(:,:,:,3).^2;%squared norm!
n2=v2(:,:,:,1).^2+v2(:,:,:,2).^2+v2(:,:,:,3).^2;%squared norm!
toc; 
clear coeff  fx fy fz


%% verifications
fprintf('Maximal norm of v1: %d\n',max(n1(:)));
fprintf('Maximal norm of v2: %d\n',max(n2(:)));
D=v1(:,:,:,1).*v2(:,:,:,1)+v1(:,:,:,2).*v2(:,:,:,2)+v1(:,:,:,3).*v2(:,:,:,3);
fprintf('Maximal scalar product between v1 and v2: %d\n', max(D(:)));
D=abs(D)>epsilon;
n2=n2~=0;
fprintf('Number of non zero v2s: %d\n', sum(n2(:)));
fprintf('Number of wrong v2s: %d\n',sum(D(:)));


mult2=(abs(dst)<epsilon).*(abs(trexp)<epsilon).*(abs(m)>epsilon)+(abs(dst)>epsilon).*(abs(p)>epsilon).*(abs(q-p.^1.5)<epsilon);
mult2=(tr>0).*mult2;
NB_mult2=sum(mult2(:));
fprintf('Cases lambda has multiplicity 2 (theory): %u\n', NB_mult2);

clear dst m trexp p q tr


if NB_mult2>0 %plotting the voxels having double multiplicity
    figure 
    v = double(f>0);
    p = patch( isosurface(v,0) ); 
    isonormals(v, p); 
    set(p, 'FaceColor','r', 'EdgeColor','none');
    axis([0 size(f,2) 0 size(f,1) 0 size(f,3)]);
    grid on;
    alpha(0.5);lighting gouraud;
    daspect([1 1 1])
    title('Points having a double orientation');
    hold on
    v = double(mult2>0);
    p = patch( isosurface(v,0) ); 
    isonormals(v, p);
    hold off
end
clear n1 n2 D trexp mult2


%% saving and 3D output
%save('test.mat');
showVectorfield3D(f,v1,v2);