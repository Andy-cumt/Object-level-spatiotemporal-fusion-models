%Guo, Dizhou 
% A Flexible Object-Level Processing Strategy to Enhance the Weight Function-Based Spatiotemporal Fusion Method
% Object-Level Hybrid Spatiotemporal Fusion: Reaching a Better Tradeoff Among Spectral Accuracy, Spatial Accuracy, and Efficiency


clear all
clc

DN_min = 0;  %% minimal DN value
DN_max = 10000;   %% maximal DN value; if using reflectance, use 1 as DN_max

Factor=16;



filename1 = 'E:\\ThirdPaper\Ukraine\L1';
FRT1 = enviread(filename1);
FRT1(FRT1<DN_min)=DN_min;

filename2 = 'E:\\ThirdPaper\Ukraine\L4';
FRT3 = enviread(filename2);
FRT3(FRT3<DN_min)=DN_min;

filename3 = 'E:\\ThirdPaper\Ukraine\M1';
CRT1 = enviread(filename3);
CRT1(CRT1<DN_min)=DN_min;

filename4 = 'E:\\ThirdPaper\Ukraine\M4';
CRT3 = enviread(filename4);
CRT3(CRT3<DN_min)=DN_min;

filename5 = 'E:\\ThirdPaper\Ukraine\M3';
CRT2 = enviread(filename5);
CRT2(CRT2<DN_min)=DN_min;


outputname='E:\已发表论文或不用的实验材料\ThirdPaper\Belle\select\M14M3OB-ESTARFM';

segfile = 'E:\已发表论文或不用的实验材料\ThirdPaper\Belle\select\segL4-L1.tif';
seg = imread(segfile);
snum = max(max(seg));

[xH,yH,bands]=size(FRT1);
xL=xH/Factor;
yL=yH/Factor;

CRT1C=zeros(xL,yL,bands);
CRT2C=zeros(xL,yL,bands);
CRT3C=zeros(xL,yL,bands);


for i=1:xL
    for j=1:yL
        tmp1 = CRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        CRT1C(i,j,:)=sum(sum(tmp1))/Factor^2;
        tmp2 = CRT2( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        CRT2C(i,j,:)=sum(sum(tmp2))/Factor^2;
        tmp3 = CRT3( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        CRT3C(i,j,:)=sum(sum(tmp3))/Factor^2;

    end
end

cubic_CRT2C = (imresize(CRT2C,Factor,'cubic'));

prediction1=zeros(xH,yH,bands);
prediction2=zeros(xH,yH,bands);
prediction=zeros(xH,yH,bands);

for b=1:bands
    f1 = FRT1(:,:,b);
    c1 = CRT1(:,:,b);
    f3 = FRT3(:,:,b);
    c3 = CRT3(:,:,b);
    c2 = CRT2(:,:,b);
    prediction1b=prediction1(:,:,b);
    prediction2b=prediction2(:,:,b);
    predictionb=prediction(:,:,b);

    for s=0:snum %注意这里从0开始
        [a1, b1, v] = find(seg == s);
        seg_pos = find(seg == s);
        num_pos = numel(v);
        c2c1=c2(seg_pos)-c1(seg_pos);
        c2c1=c2c1(:);
        mean_c21=mean(c2c1);
        c2c3=c2(seg_pos)-c3(seg_pos);
        c2c3=c2c3(:);
        mean_c23=mean(c2c3);
        
        c1c3=c1(seg_pos)-c3(seg_pos);
        c1c3=c1c3(:);
        mean_c13=abs(mean(c1c3));
        
        T_1=1/abs(sum(c2c1)+0.0001);
        T_2=1/abs(sum(c2c3)+0.0001);
        T1=T_1/(T_1+T_2);
        T2=T_2/(T_1+T_2);
        
        x1=f1(seg_pos);
        x3=f3(seg_pos);
        x=[x1;x3];
        
        y1=c1(seg_pos);
        y3=c3(seg_pos);
        y=[y1;y3];
        
        xx=[ones(length(y),1),x];
        [BB,~,~,~,stats]=regress(y,xx);
        Gain=BB(2);
%        bias=BB(1);
        clear x1 y1 x3 y3 x y xx
%        test(s+1)=Gain;
        %如果可供回归的像元很少，也要考虑
        
        if (Gain>0) && (Gain<5) && (mean_c13 > 0.02*DN_max) && (stats(3) <= 0.05)
            Gain=BB(2);
        else
            Gain=1.0;
        end
        

        prediction1b(seg_pos)=f1(seg_pos)+Gain*mean_c21;
        prediction2b(seg_pos)=f3(seg_pos)+Gain*mean_c23;
        predictionb(seg_pos)=T1*prediction1b(seg_pos)+T2*prediction2b(seg_pos);
        
    end
    prediction(:,:,b)=predictionb(:,:);
end
% predictionC1=zeros(xL,yL,bands);
% predictionC2=zeros(xL,yL,bands);
% FRT1C=zeros(xL,yL,bands);
% FRT3C=zeros(xL,yL,bands);
% for i=1:xL
%     for j=1:yL
%         tmp1 = prediction1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
%         predictionC1(i,j,:)=sum(sum(tmp1))/Factor^2;
%         tmp2 = prediction2( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
%         predictionC2(i,j,:)=sum(sum(tmp2))/Factor^2;
%         tmp3 = FRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
%         FRT1C(i,j,:)=sum(sum(tmp3))/Factor^2;
%         tmp4 = FRT3( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
%         FRT3C(i,j,:)=sum(sum(tmp4))/Factor^2;
%     end
% end
% residual1=(CRT2C-CRT1C)-(predictionC1-FRT1C);
% cubic_residual1=(imresize(residual1,Factor,'cubic'));
% residual2=(CRT2C-CRT3C)-(predictionC2-FRT3C);
% cubic_residual2=(imresize(residual2,Factor,'cubic'));
% prediction=prediction+T1*cubic_residual1+T2*cubic_residual2;



    for b=1:bands
        for i=1:xH
            for j=1:yH
                if prediction(i,j,b)<DN_min
                   prediction_pos=max(DN_min,cubic_CRT2C(i,j,b));
                   prediction(i,j,b)=min(DN_max,prediction_pos);
                end
                if prediction(i,j,b)>DN_max
                   prediction_pos=min(DN_max,cubic_CRT2C(i,j,b));
                   prediction(i,j,b)=max(DN_min,prediction_pos);
                end
            end
        end
    end



    
file_patch = prediction;
info=enviinfo(file_patch);
enviwrite(file_patch,info,outputname);  











