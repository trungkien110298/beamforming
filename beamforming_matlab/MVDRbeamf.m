function y = MVDRbeamf(FX,PhiNN,X,H,post_filtering)
% input
%FX = num_mic x nFreq
%PhiNN = num_mic x num_mic x nFreq
%X = num_mic X nFrame x nFreq
[num_mic,nFrame,nFreq] = size(X);
channel_selection = zeros(num_mic,1);
smooth_post_filter = zeros(nFreq,nFrame);

for i = 1:nFrame
    for j = 1: nFreq
        r_inv = inv(PhiNN(:, :,j) + (1e-4) * diag(ones(num_mic, 1))); 
        
        w(:, i,j) = r_inv * FX(:,j) ./ (conj(FX(:,j).') * r_inv * FX(:,j)); % MVDR        
        y(i,j) = (w(:, i,j))'*X(:,i,j);
%         % post filtering
%         phiss(i,j) = y(i,j)*y(i,j);
%         phinn(i,j) = 1./((conj(FX(:,j).') * r_inv * FX(:,j)));
%         post_filter(i,j) = phiss(i,j)./(phiss(i,j)+phinn(i,j));
%         mask(i,j) = min(abs(post_filter(i,j)),1);
    end
%     smooth_post_filter(1:nFreq-1,i) = gammatone_smooth(mask(i,:),abs(y(i,:)),H);
%     smooth_post_filter(nFreq,i) = 1e-4;
%     if strcmp(post_filtering,'yes')
%         y(i,:) = smooth_post_filter(:,i)'.*y(i,:);
%     else
%         y(i,:) = y(i,:);
%     end
end
end





