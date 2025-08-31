%Guo, Dizhou 
% A Flexible Object-Level Processing Strategy to Enhance the Weight Function-Based Spatiotemporal Fusion Method
% Object-Level Hybrid Spatiotemporal Fusion: Reaching a Better Tradeoff Among Spectral Accuracy, Spatial Accuracy, and Efficiency



clear all
clc

DN_min = 0;  %% minimal DN value
DN_max = 10000;   %% maximal DN value; if using reflectance, use 1 as DN_max

Factor=30;

%similar_pixel=50;

%%%%%%%%%

filename1 = 'C:\Users\Dizhou Guo\Downloads\S2_2021-06-11_SR_10000.tif';
FRT1 = enviread(filename1);
FRT1(FRT1<DN_min)=DN_min;

filename2 = 'C:\Users\Dizhou Guo\Downloads\S3_2021-06-11_SR_Final_10000.tif';
CRT1 = enviread(filename2);
CRT1(CRT1<DN_min)=DN_min;

filename3 = 'C:\Users\Dizhou Guo\Downloads\S3_2021-07-16_SR_Final_10000.tif';
CRT2 = enviread(filename3);
CRT2(CRT2<DN_min)=DN_min;

outputname='C:\Users\Dizhou Guo\Downloads\S2S3';

segfile = 'C:\Users\Dizhou Guo\Downloads\Segmentation_Result.tif';
seg = imread(segfile);
snum = max(max(seg));


CRT2(CRT2<DN_min)=DN_min;

[xH,yH,bands]=size(FRT1);

CRT1=double(CRT1);
CRT2=double(CRT2);
FRT1=double(FRT1);

xL=xH/Factor;
yL=yH/Factor;

CRT1C=zeros(xL,yL,bands);
CRT2C=zeros(xL,yL,bands);

for i=1:xL
    for j=1:yL
        tmp1 = CRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        CRT1C(i,j,:)=sum(sum(tmp1))/Factor^2;
        tmp2 = CRT2( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        CRT2C(i,j,:)=sum(sum(tmp2))/Factor^2;

    end
end


cubic_CR = (imresize(CRT2C-CRT1C,Factor,'cubic'));


%%%%%%%%%%%%%%%%   downscale  %%%%%%%%%%%%%%%%%

prediction=zeros(xH,yH,bands);





for b=1:bands
    f1 = FRT1(:,:,b);
    c1 = CRT1(:,:,b);        
    c2 = CRT2(:,:,b);
    predictionb=prediction(:,:,b);

    for s=0:snum %注意这里从0开始
        [a1, b1, v] = find(seg == s);
        seg_pos = find(seg == s);
        num_pos = numel(v);
        c2c1=c2(seg_pos)-c1(seg_pos);
        c2c1=c2c1(:);
        median_c=median(c2c1);
        predictionb(seg_pos)=median_c+f1(seg_pos);    

    end
    prediction(:,:,b)=predictionb(:,:);
end
predictionC=zeros(xL,yL,bands);
FRT1C=zeros(xL,yL,bands);
for i=1:xL
    for j=1:yL
        tmp1 = prediction( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        predictionC(i,j,:)=sum(sum(tmp1))/Factor^2;
        tmp2 = FRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        FRT1C(i,j,:)=sum(sum(tmp2))/Factor^2;
    end
end
residual=(CRT2C-CRT1C)-(predictionC-FRT1C);
cubic_residual=(imresize(residual,Factor,'cubic'));
prediction=prediction+cubic_residual;



    for b=1:bands
        for i=1:xH
            for j=1:yH
                if prediction(i,j,b)<DN_min
                   prediction_pos=max(DN_min,FRT1(i,j,b)+cubic_CR(i,j,b));
                   prediction(i,j,b)=min(DN_max,prediction_pos);
                end
                if prediction(i,j,b)>DN_max
                   prediction_pos=min(DN_max,FRT1(i,j,b)+cubic_CR(i,j,b));
                   prediction(i,j,b)=max(DN_min,prediction_pos);
                end
            end
        end
    end


file_patch = prediction;
info=enviinfo(file_patch);
enviwrite(file_patch,info,outputname);  


