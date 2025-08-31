function [a,b] = guidedfilter_MS_low_Dizhou(I, p)

%   GUIDEDFILTER_MS   O(1) time implementation of guided filter using a MS image as the guidance.
%
%   - guidance image: I (should be a MS image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps

[hei, wid] = size(p);
dim=size(I,3);
N = hei*wid; % the size of each local patch; N=(2r+1)^2 except for boundary pixels.这里是提高窗口大小

mean_p = mean2(p);%每个窗口中的M2平均值
mean_I= mean2(I);%%每个窗口中的M1平均值
mean_Ip=mean2(I.*p);
cov_Ip= mean_Ip- mean_I*mean_p;% covariance of (I, p) in each local patch.


% variance of I in each local patch: the matrix Sigma in Eqn (14).
% Note the variance in each local patch is a dim x dim symmetric matrix:

        var_I=mean2(I.*I)  - mean_I*mean_I;  
 
%a = zeros(hei, wid, dim);
%for y=1:hei
%    for x=1:wid
        %for k=1:dim  Sigma(k,:)=var_I(y,x,k,:);  end
        Sigma=reshape(var_I,dim,dim);
        %cov_Ip0=D3_D2(cov_Ip)';        
        a = cov_Ip * inv(Sigma + eps * eye(dim)); % Eqn. (14) in the paper;
%    end
%end
b = mean_p - sum(a.* mean_I,3); % Eqn. (15) in the paper;

% q=boxfilter(b, r);
% for i=1:dim
%     q=q+boxfilter(a(:,:,i), r).* I(:, :,i);
% end
% q=q./N;
end