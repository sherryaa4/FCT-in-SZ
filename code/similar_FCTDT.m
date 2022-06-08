function [] = similar_FCTDT(FCT_path,DTI_path,mask_path,subID)
% Summary of this function goes here 
% computing similarity between Functional Correlation Tensor and diffusion tensor

% Input
% FCT_path:          the functional correlation tensor path 
% DTI_path:          the DTI data path (must be calculated by MRtrix3)
% mask_path:         the mask image of FCT/DT image
% subID:             the subject ID 


% Note!!!!!: the functional correlation tensor and diffusion tensor must be saved in same format (fsl or mrtrix3) and in same space (diffusion space or MNI space)

% Written by Jiajia Zhao
% /2022/06/08

% the equation used to calculated the similarity between functional correlation tensor and diffusion tensor was revised by following referrence:
% Alexander, D., Gee, J., & Bajcsy, R. (1999). Similarity Measures for Matching Diffusion Tensor Images Procedings of the British Machine Vision Conference 1999.




%--------------------------------Main Function---------------------------------------%



    
        fprintf('Processing  : Computing Similarity between FCT and DT...\n.')
        T1=load_untouch_nii(FCT_path);
        T2=load_untouch_nii(DTI_path);
        
        brain_mask=load_untouch_nii(mask_path);
        Img_mask=brain_mask.img;
        ind=find(Img_mask~=0 & ~isnan(Img_mask));
        Mask_xyz=zeros(length(ind),3);
        [Mask_xyz(:,1),Mask_xyz(:,2),Mask_xyz(:,3)]=ind2sub(size(Img_mask),ind);
        nvox=size(Mask_xyz,1);
        
        [~,d1]=tensor_vec2mat_mrtrix(T1.img,Img_mask);
        [~,d2]=tensor_vec2mat_mrtrix(T2.img,Img_mask);
        
        simi_d=NaN(nvox,1);
        parfor n=1:nvox
            
            if sum(sum(d1{n}==0))~=9 && sum(sum(d2{n}==0))~=9 && ~isempty(d1{n}) && ~isempty(d2{n})
                dif=d1{n}-d2{n};
                euc_dis=sqrt(trace(dif^2));
                simi_d(n,1)=1/euc_dis;
            end
            
        end
        indinf=simi_d==Inf;
        simi_d(indinf)=nan;
        
        B2=NaN(size(Img_mask));
        for n=1:nvox
            B2(Mask_xyz(n,1),Mask_xyz(n,2),Mask_xyz(n,3))=simi_d(n);
        end
        brain_mask.img=B2;
        brain_mask.hdr.dime.datatype=16;
        brain_mask.hdr.dime.bitpix=32;
        filename=[subID,'_DTIFCT_similar.nii.gz'];
        save_untouch_nii(brain_mask,filename)
        
    end
end
