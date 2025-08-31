%Guo, Dizhou 
% A Flexible Object-Level Processing Strategy to Enhance the Weight Function-Based Spatiotemporal Fusion Method
% Object-Level Hybrid Spatiotemporal Fusion: Reaching a Better Tradeoff Among Spectral Accuracy, Spatial Accuracy, and Efficiency

clear all
clc

DN_min = 0;  %% minimal DN value
DN_max = 10000;   %% maximal DN value; if using reflectance, use 1 as DN_max

Factor=8;
Bands = 6;

tif='.tif';
tiff='.tiff';


path = 'E:\deeplearningdata\newcolby\test\20190217' ;
filename1 = fullfile(path,'0_LC08_031033_20190217cutt.tif');
filename2 = fullfile(path,'0_MCD_20190217cutt.tif');
filename3 = fullfile(path,'1_MCD_20190828cutt.tif');
segfile = fullfile(path,'seg_LC08_031033_20190217cutt.tif');
outname=fullfile(path,'OLHSTFM_20190828');



[~,~,format]=fileparts(filename1);
if strcmp(format,tif) || strcmp(format,tiff)
    FRT1 = imread(filename1);
else
    FRT1 = enviread(filename1);
end
FRT1(FRT1<DN_min)=DN_min;
FRT1(FRT1>DN_max)=DN_max;


[~,~,format]=fileparts(filename2);
if strcmp(format,tif) || strcmp(format,tiff)
    CRT1 = imread(filename2);
else
    CRT1 = enviread(filename2);
end
CRT1(CRT1<DN_min)=DN_min;
CRT1(CRT1>DN_max)=DN_max;
CRT1 = (imresize(CRT1,Factor,'nearest'));


[~,~,format]=fileparts(filename3);
if strcmp(format,tif) || strcmp(format,tiff)
    CRT2 = imread(filename3);
else
    CRT2 = enviread(filename3);
end
CRT2(CRT2<DN_min)=DN_min;
CRT2(CRT2>DN_max)=DN_max;
CRT2 = (imresize(CRT2,Factor,'nearest'));


FRT1 = FRT1(:,:,1:Bands);
CRT1 = CRT1(:,:,1:Bands);
CRT2 = CRT2(:,:,1:Bands);


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

cubic_CRT1 = (imresize(CRT1C,Factor,'cubic'));
cubic_CR = (imresize(CRT2C-CRT1C,Factor,'cubic'));

edge = zeros(xH,yH,bands);
for b=1:bands
    edge(:,:,b)=Robert(FRT1(:,:,b),0.3);
end

%%%%%%%%%%%%%%%%   downscale  %%%%%%%%%%%%%%%%%

prediction_F=zeros(xH,yH,bands);
prediction_S=zeros(xH,yH,bands);

mean_FRT1_allb=mean(mean(FRT1,1),2);
mean_FRT1_allb(:)=mean_FRT1_allb;
FRT1b=zeros(xH,yH,bands);

RI=zeros(xH,yH,bands)+1;


%计算中心和f1平均值
central_seg_x = zeros(snum+1,1);
central_seg_y = zeros(snum+1,1);
f1_seg_mean = zeros(snum+1,bands);
num_pos_seg = zeros(snum+1,1);

for s=0:snum %注意这里从0开始
    [a1, b1, v] = find(seg == s);
    seg_pos = find(seg == s);
    num_pos_seg(s+1) = numel(v);
    central_seg_x(s+1)=mean(a1);
    central_seg_y(s+1)=mean(b1);
    for b=1:bands
        f1 = FRT1(:,:,b);
        RIb=RI(:,:,b);
        f1_seg_mean(s+1,b)=mean(f1(seg_pos));
       
    end
end
Gain_image=zeros(xH,yH,bands);
bias_image=zeros(xH,yH,bands);
abnumber_image=zeros(xH,yH,bands);
guidedfilterA = zeros(xH,yH,bands);
OLedge = zeros(xH,yH,bands);
OLedgenew = zeros(xH,yH,bands);
for b=1:bands
    f1 = FRT1(:,:,b);
    c1 = CRT1(:,:,b);
    c2 = CRT2(:,:,b);
    cc1=cubic_CRT1(:,:,b);
    cc=cubic_CR(:,:,b);
    bias_imageb=bias_image(:,:,b);
    Gain_imageb=Gain_image(:,:,b);
    f1_seg_meanb=f1_seg_mean(:,b);
    predictionFb=prediction_F(:,:,b);
    predictionSb=prediction_S(:,:,b);
    guidedfilterAb = guidedfilterA(:,:,b);
    edgeb = edge(:,:,b);
    
    
    edgenew = abs(edgeb-mean2(edgeb))./(4.*std2(edgeb));
    edgenew(edgenew>1) = 1;
    
    for s=0:snum %注意这里从0开始
        [a1, b1, v] = find(seg == s);
        seg_pos = find(seg == s);
        num_pos = numel(v);
        c2c1=c2(seg_pos)-c1(seg_pos);
        c2c1=c2c1(:);
        num_diff = length(unique(c2c1));
        median_c=median(c2c1);
        predictionSb(seg_pos)=median_c+f1(seg_pos);   
      
        
        if (num_diff < 3)||(num_pos < 4*Factor^2)
            central_seg_x0 = central_seg_x(s+1); 
            central_seg_y0 = central_seg_y(s+1);
            
            w=3*Factor;
            ai=max(1,round(central_seg_x0-w));
            bi=min(xH,round(central_seg_x0+w));
            aj=max(1,round(central_seg_y0-w));
            bj=min(yH,round(central_seg_y0+w));
            %用能把块包裹住的最小四边形
            ai2=min(a1);
            bi2=max(a1);
            aj2=min(b1);
            bj2=max(b1);
            ai=min(ai,ai2);
            bi=max(bi,bi2);
            aj=min(aj,aj2);
            bj=max(bj,bj2);
            w_image=zeros(xH,yH);
            w_image(ai:bi,aj:bj)=1;
            w_index=find(w_image == 1);
                    
            similar_index=w_index;

            num_similar=length(similar_index(:));
            num=num_similar;

                seg_pos3=union(similar_index,seg_pos);
                
                x=c1(ai:bi,aj:bj);
                y=c2(ai:bi,aj:bj);
                
                [Gain,bias] = guidedfilter_MS_low_Dizhou(x, y);

                  clear x y
            predictionFb(seg_pos)=Gain*f1(seg_pos)+bias;
            guidedfilterAb(seg_pos)=Gain;
               clear seg_pos3
        else
                x(:)=c1(seg_pos);
                y(:)=c2(seg_pos);

                [Gain,bias] = guidedfilter_MS_low_Dizhou(x, y);
                 clear x y 
                predictionFb(seg_pos)=Gain*f1(seg_pos)+bias;
                guidedfilterAb(seg_pos)=Gain;

        end
    end
     prediction_F(:,:,b)=predictionFb(:,:);  
     prediction_S(:,:,b)=predictionSb(:,:);
     guidedfilterA(:,:,b)=guidedfilterAb(:,:);
    
    OLedgenew(:,:,b) = edgenew(:,:);
end

% 
predictionFC=zeros(xL,yL,bands);
predictionSC=zeros(xL,yL,bands);
FRT1C=zeros(xL,yL,bands);
for i=1:xL
    for j=1:yL
        tmp1 = prediction_F( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        predictionFC(i,j,:)=sum(sum(tmp1))/Factor^2; 
        tmp2 = prediction_S( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        predictionSC(i,j,:)=sum(sum(tmp2))/Factor^2; 
        tmp3 = FRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        FRT1C(i,j,:)=sum(sum(tmp3))/Factor^2;
    end
end
clear tmp1 tmp2 tmp3
nearest_prediction=(imresize(predictionFC,Factor,'nearest'));

residual=(CRT2C-CRT1C)-(predictionSC-FRT1C);
cubic_residual=(imresize(residual,Factor,'cubic'));
prediction_S=prediction_S+cubic_residual;

for b=1:bands
    nearest_predictionb=nearest_prediction(:,:,b);
    c2 = CRT2(:,:,b);    
    for s=0:snum %注意这里从0开始
        [a1, b1, v] = find(seg == s);
        seg_pos = find(seg == s);
        num_pos = numel(v);
        change=mean2(c2(seg_pos)-nearest_predictionb(seg_pos));
        for ii = 1:num_pos
            a0 = a1(ii);
            b0 = b1(ii);                        
            prediction_F(a0,b0,b)=change+prediction_F(a0,b0,b);
        end
    end
end



OLedgenew = (imresize(OLedgenew,0.5,'nearest'));
OLedgenew = (imresize(OLedgenew,2,'cubic'));

T = 0.1;
weight = zeros(xH,yH,bands);
FitFCA = abs(guidedfilterA);

index = find(FitFCA<1 & OLedgenew>T);



prediction = prediction_F;
weight(index) = (1-FitFCA(index)).*(1-OLedgenew(index)+T);
prediction(index) = prediction_S(index).* weight(index) + prediction_F(index).*(1-weight(index));

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
enviwrite(file_patch,info,outname);  