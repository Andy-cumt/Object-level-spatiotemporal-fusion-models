function OutImg = Normalize(InImg)
    ymax=255;ymin=0;
    xmax = max(max(InImg)); %���InImg�е����ֵ
    xmin = min(min(InImg)); %���InImg�е���Сֵ
    OutImg = round((ymax-ymin)*(InImg-xmin)/(xmax-xmin) + ymin); %��һ����ȡ��
end
