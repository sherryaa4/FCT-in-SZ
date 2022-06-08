function [neib_vox_vec,vox_neib_xyz,vox_xyz] = neib_vec2(brain_mask)
%neib_vox Summary of this function goes here
% find 26-neibours of each voxel

%   Detailed explanation goes here

%  INPUT
%  brain_mask:         a 3d .img

%  OUTPUT

%  neib_vox_vec:      the vector between each voxel and its neibours, a cell matrix with total number of 
%                     voxels(n_x*n_y*n_z) x 1 dimension, in the cell contains 3d unit vector between each voxel 
%                     and its neibours, a matrix with number of neibours x 3 dimension

%  neib_vox_ind:      the index of each voxel's neibours, a matrix with total number of voxels(n_x*n_y*n_z) x 26 dimension, the coordinate 
%                     of voxel was transfered into index coordinate, not 3d coordinate
 
%----------------------------------------------------------------------------------------------------------------------------------%


[vox_xyz(:,1),vox_xyz(:,2),vox_xyz(:,3)]=ind2sub(size(brain_mask),find(brain_mask~=0));

vox_x=vox_xyz(:,1);
vox_y=vox_xyz(:,2);
vox_z=vox_xyz(:,3);
num_vox=size(vox_xyz,1);
vox_neib_xyz{num_vox,1}=[];
neib_vox_vec{num_vox,1}=[];
% the 26-neibours

parfor v=1:num_vox
      
%           the 26 neibours
            t=1;
            tmp_neib_vox=zeros(27,3);
            tmp_neib_vec=zeros(27,3);
            for a=-1:1   
           %for a=1:-1:-1    
                for b=-1:1
                    for c=-1:1
                        tmp_neib_vox(t,:)=[vox_x(v)-a,vox_y(v)+b,vox_z(v)+c];
                        tmp_sqr=sqrt(a*a+b*b+c*c);
                        tmp_neib_vec(t,:)=[a/tmp_sqr,b/tmp_sqr,c/tmp_sqr];
                        
                        t=t+1;
                    end
                end
            end
            tmp_neib_vox(14,:)=[]; % remove itself while a=b=c=0
            tmp_neib_vec(14,:)=[];
            
            vox_neib_xyz{v,1}=tmp_neib_vox;
            neib_vox_vec{v,1}=tmp_neib_vec;

            
end


end