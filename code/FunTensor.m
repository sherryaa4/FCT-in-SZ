function [] = FunTensor(fmri_path,mask_path,subID)
% Summary of this function goes here 
% computing Functional Correlation Tensor

% Input
% fmri_path:         the fmri data path
% mask_path:         the mask file for the fmri data
% subID:                the subject ID 

% Output


% Written by Jiajia Zhao
% /2022/06/08

% The detailed interpretation of method to calculate the functional tensor can be found in references as follows 

% Ding, Z., Newton, A. T., Xu, R., Anderson, A. W., Morgan, V. L., & Gore, J. C. (2013). Spatio-temporal correlation tensors reveal functional structure in human brain. PloS one, 8(12), e82107.
% Ding, Z., Xu, R., Bailey, S. K., Wu, T. L., Morgan, V. L., Cutting, L. E., ... & Gore, J. C. (2016). Visualizing functional pathways in the human brain using correlation tensors and magnetic resonance imaging. Magnetic resonance imaging, 34(1), 8-17.





%--------------------------------Main Function---------------------------------------%




mask=load_untouch_nii(mask_path);
brain_mask=mask.img;
[vox_dimenx,vox_dimeny,vox_dimenz]=size(brain_mask);
[neib_vox_vec,vox_neib_xyz,vox_xyz] = neib_vec2(brain_mask);
num_vox=size(vox_xyz,1);

    fprintf('Processing  : Functional Correlation Tensor Computing....')
    sub_nii=load_untouch_nii(fmri_path);
    rest=sub_nii.img;
    
    [neib_cor] = NeibCor2(vox_xyz,vox_neib_xyz,rest);
    zneib_cor=fisherz(neib_cor);
    C=zneib_cor.*zneib_cor;
    T=zeros(num_vox,6);
   
  
    parfor i=1:num_vox
        tmp_vec=neib_vox_vec;
        neib_vec=tmp_vec{i,1};
        M=designM(neib_vec);
        M_transp=M';
        neib_C=C(i,:)';
        tmp_T=(inv(M_transp*M))*M_transp*neib_C;
        T(i,:)=tmp_T';
    end
    

%% create tensor matrix for MRtrix

% designM/T: xx, 2xy, 2xz, yy, 2yz, zz


%The tensor coefficients are stored in the output image as follows:
%volumes 0-5: D11, D22, D33, D12, D13, D23  MRview

vol1=T(:,1);
vol2=T(:,4);
vol3=T(:,6);
vol4=T(:,2);
vol5=T(:,3);
vol6=T(:,5);



B=[vol1,vol2,vol3,vol4,vol5,vol6];
B2=zeros(vox_dimenx,vox_dimeny,vox_dimenz,6);
for n=1:num_vox
    B2(vox_xyz(n,1),vox_xyz(n,2),vox_xyz(n,3),:)=B(n,:);
end
sub_nii.img=B2;
sub_nii.hdr.dime.dim(5)=6;
sub_nii.hdr.dime.pixdim(5)=1;
filename=[subID,'_inDWIspace_FCT.nii.gz'];
save_untouch_nii(sub_nii,filename)



end