function [M] = designM(neib_vec)
% dyadic_tensor Summary of this function goes here
% calculate design matrix
%   Detailed explanation goes here

%INPUT
% neib_vec:    unit vector for 26 neighbors, a matrix with 26 X 3 dimension

%OUTPUT
% M:           neibour design matrix , a matrix with 26 X 6 dimension


[Nneib,~]=size(neib_vec);
M=zeros(Nneib,6);
for neib=1:Nneib
  M(neib,:)=dyadic_tensor_half(neib_vec(neib,:));
end


end


function [D_tensor] = dyadic_tensor_half(unit_vec)
% dyadic_tensor Summary of this function goes here
% calculate dyadic tensor
%   Detailed explanation goes here

%INPUT
% unit_vec:    a 3d unit vector, a matrix with 1 X 3 dimension

%OUTPUT
% D_tensor:    dyadic tensor, a matrix with 1 X 6 dimension:(x*x,2*x*y,2*x*z,y*y,2*y*z,z*z)
% 

vec_x=unit_vec(1);
vec_y=unit_vec(2);
vec_z=unit_vec(3);
D_tensor=zeros(1,6);
D_tensor(1,1)=vec_x*vec_x;
D_tensor(1,2)=2*vec_x*vec_y;
D_tensor(1,3)=2*vec_x*vec_z;
D_tensor(1,4)=vec_y*vec_y;
D_tensor(1,5)=2*vec_y*vec_z;
D_tensor(1,6)=vec_z*vec_z;
% D_tensor(1,1)=vec_x*vec_x;
% D_tensor(1,2)=vec_x*vec_y;
% D_tensor(1,3)=vec_x*vec_z;
% D_tensor(1,4)=vec_y*vec_y;
% D_tensor(1,5)=vec_y*vec_z;
% D_tensor(1,6)=vec_z*vec_z;


end