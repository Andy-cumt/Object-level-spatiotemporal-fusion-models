clear all
clc

DN_min = 0;  %% minimal DN value
DN_max = 10000;   %% maximal DN value; if using reflectance, use 1 as DN_max

Factor=20;

%similar_patch=20;

%%%%%%%%%

filename1 = 'E:\文档迁移\RASDF融合\异质型800最后论文2\无误差\真实影像\新加实验\M1M5OBLRSFMab';
prediction = enviread(filename1);
prediction(prediction<DN_min)=DN_min;

filename3 = 'E:\文档迁移\RASDF融合\异质型800最后论文2\无误差\真实影像\新加实验\M5';
CRT2 = enviread(filename3);
CRT2(CRT2<DN_min)=DN_min;

segfile = 'E:\文档迁移\RASDF融合\异质型800最后论文2\无误差\重采样影像\seg.tif';
seg = imread(segfile);
snum = max(max(seg));

[xH,yH,bands]=size(prediction);
xL=xH/Factor;
yL=yH/Factor;
predictionC=zeros(xL,yL,bands);

for i=1:xL
    for j=1:yL
        tmp1 = prediction( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        predictionC(i,j,:)=sum(sum(tmp1))/Factor^2;        
    end
end
nearest_prediction=(imresize(predictionC,Factor,'nearest'));

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
            prediction(a0,b0,b)=change+prediction(a0,b0,b);
        end
    end
end
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
enviwrite(file_patch,info,'E:\文档迁移\RASDF融合\异质型800最后论文2\无误差\真实影像\新加实验\M1M5OBLRSFM');  