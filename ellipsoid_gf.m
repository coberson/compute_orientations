function I=ellipsoid_gf(c_x,c_y,c_z,l_1,l_2,l_3,a_1,a_2,n_x,n_y,n_z)
% gives the 3D matrix of an ellipsoid (general form) image:
% c_x,c_y, c_z is the center
% l_1, l_2, l_3 are the semi-axes
% a_1 and a_2 are the angles defining v_1, v_2 and v_3 the principal axes
% of the ellipse
% n_x,n_y,n_z are the half dimensions of the image
% I is the ellipsoid image
v_1=[cos(a_1)*sin(a_2);sin(a_1)*sin(a_2);cos(a_2)];
v_2=[cos(a_1)*cos(a_2);sin(a_1)*cos(a_2);-sin(a_2)];
v_3=cross(v_1,v_2);
[x,y,z] = ndgrid(-n_x:n_x,-n_y:n_y,-n_z:n_z);
A=diag([1/l_1^2,1/l_2^2,1/l_3^2]);
R=[v_1,v_2,v_3];
K=R*A*R';
I=(x-c_x).^2*K(1,1)+(x-c_x).*(y-c_y)*K(2,1)+(x-c_x).*(z-c_z)*K(3,1)+(y-c_y).*(x-c_x)*K(1,2)+(y-c_y).^2*K(2,2)+(y-c_y).*(z-c_z)*K(3,2)+(z-c_z).*(x-c_x)*K(1,3)+(z-c_z).*(y-c_y)*K(2,3)+(z-c_z).^2*K(3,3)<=1;



% figure
% v = double(I>0);
% p = patch( isosurface(v,0) );
% isonormals(v, p);
% set(p, 'FaceColor','r', 'EdgeColor','none');
% axis([0 size(I,2) 0 size(I,1) 0 size(I,3)]);
% grid on;
% alpha(0.5);lighting gouraud;
% daspect([1 1 1])
% title('Ellipsoid');

