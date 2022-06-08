function [neib_cor] = NeibCor2(vox_xyz,vox_neib_xyz,rest)
%NeibCor Summary of this function goes here
% calculate correlation between each voxel and its neibours
%   Detailed explanation goes here

%INPUT 
%  vox_xyz:           the x,y,z coordinate of each voxel

%  vox_neib_xyz:      the x,y,z coordinate of each voxel's neibours, a cell matrix with total number of voxels
%                     (n_x*n_y*n_z) x 1 dimension, in the cell contains x,y,z coordinate of each voxel's neibours, 
%                      a matrix with number of neibours x 3 dimension  
%  rest:              the 4d image, a matrix with n_x*n_y*n_z*n_tc dimension,including total dimension of x,y,z axis and length of timeseries 

%OUTPUT
%  neib_cor:          the correlation between each voxel and its neibours, a matrix with total number of voxels(n_x*n_y*n_z) x 26 dimension

%----------------------------------------------------------------------------------------------------------------------------------%

vox_x=vox_xyz(:,1);
vox_y=vox_xyz(:,2);
vox_z=vox_xyz(:,3);
[nx,ny,nz,n_len]=size(rest);
%nz=size(rest);
n_vox=size(vox_xyz,1);
neib_cor=zeros(n_vox,26); 
parfor i=1:n_vox
    tmp_rest=rest;
    tc=squeeze(tmp_rest(vox_x(i),vox_y(i),vox_z(i),:));
    if sum(tc==0)==n_len     % skip the zero timeseries
        continue;
    end
    for j=1:26
        neib_xyz=vox_neib_xyz{i,1}(j,:);
        %if isempty(find(neib_xyz<1|neib_xyz>nz,3))   % skip the unlogical neib index which is out  of  dimension
        if neib_xyz(1)>0 && neib_xyz(1) <nx+1 && neib_xyz(2)>0 && neib_xyz(2) <ny+1 && neib_xyz(3)>0 && neib_xyz(3) <nz+1
            neib_tc=squeeze(tmp_rest(neib_xyz(1),neib_xyz(2),neib_xyz(3),:));
            if sum(neib_tc==0)~=n_len  % skip the zero timeseries
                tmp=corr(tc,neib_tc);
                neib_cor(i,j)=tmp;
            end
        end
        
    end
end


end
