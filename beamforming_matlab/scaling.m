function [A, YA] = scaling(X, Y, param)

%   X: sensor observations, nMic * nFrame * nFreq
%   Y: separated signals, nOut * nFrame * nFreq
%   A: output, nMic * nOut * nFreq 
%   YA: output, nMic * nFrame * nFreq * nOut (if needed)

[nMic, nFrame, nFreq] = size(X);
[nOut, nFrame, nFreq] = size(Y);

A = zeros(nMic, nOut, nFreq);
if nargout >= 2
  YA = zeros(nMic, nFrame, nFreq, nOut);
end


  
for f=freqBinRange(param)
  x = X(:,:,f);
  y = Y(:,:,f);

  A(:,:,f) = (x*y')/(y*y');
  if nargout >=2
    for i=1:nOut
      YA(:,:,f,i) = A(:,i,f)*y(i,:); 
    end
  end
end
