function speech_enhancement_derev_offline_viet

% input
% input: 8 channel audio data
%------------------------------------------------------
% preparing for copying from server to the local machine
% then process and copy it back to the sever
% addpath('wpe_v1.31/');
addpath('/home/kienpt/Documents/Beam/data/');
%scp = ['scp '];
%server_noisy_data = ['tdung@sp49.rlab.miniy.yahoo.co.jp:'];
%processed_data_dir = ['/home/tdung/data/processed/'];
%current_dir = ' .';
% next file where the wavfile is located in the wavlist file
%file_list = 'wavlist.cmd.exRS';
%fid = fopen(file_list,'rt');
%list_file = textscan(fid,'%s', 'Delimiter', '\n');

safe_region = 10;
subnoisecov = 'yes';
derev        = 'yes';

% iNumFilts = 40;
% nFFT = 512;
% H = ComputeFilterResponse(iNumFilts, nFFT);

    filename = '/home/kienpt/Documents/Beam/data/speaker0250-0020.wav';
    [x,fs] = audioread(filename);
    tic
    noise_begin = 100;
    noise_end = 0;
    for j = 1:2 
        s(:,j) = x(:,j);
        % apply VAD to find the speech boundary
        result(:,j) =  VAD(s(:,j),fs, 0.1, 0.05, 0.025, 20, 1);
        left_side = (min(find(result(:,j)==1))-safe_region)*0.025;
        left_side = max(left_side,0);
        
        right_side = (max(find(result(:,j)==1))+safe_region)*0.025;
        right_side = min(right_side,length(s(:,j))*(1/fs));
        
        % found the smallest and largest index and pick them. That is safe
        % to avoid learn noise covariance matrix in the speech region
        
        noise_begin = min(noise_begin, left_side);
        noise_end = max(noise_end,right_side);
    end
    t1 = toc;
    fprintf('VAD duration: %f\n',t1);
    % remove all wav file
    if isempty(noise_begin)
        noise_begin = 0.1;
    end
    if isempty(noise_end)
        noise_end = length(s(:,j))*(1/fs) - 0.1;
    end
    %cmd = ['rm *.wav'];
    %system(cmd);
    % note that take the last part
    noise_end = length(s(:,j))*(1/fs) - noise_end;
    % apply speech enhancement
    [mvdr_out,mvdr_derev_out] = SE(s, fs, noise_begin, noise_end,subnoisecov);
    % save data  
    audiowrite('/home/kienpt/Documents/Beam/data/speaker0250-0020_BF.wav',mvdr_out,fs);
    clear s
    clear result
    clear x
    t2 = toc;
    fprintf('SE duration: %f\n',t2 - t1);
end


function [mvdr_out,mvdr_derev_out] = SE(s,fs, NoiseBeg, NoiseEnd,subnoisecov)

tic
mvdr_out = zeros(size(s(:,1)));
mvdr_derev_out = zeros(size(s(:,1)));
% [s,fs] = audioread(file_name);
param.rate = fs;
param.freqRange = [0 fs/2];
param.tNoiseBeg = NoiseBeg;
param.tNoiseEnd = NoiseEnd;

% perform dereverberation in advance
[nb_samples,nb_chan] = size(s);
% for i = 1:nb_chan
%     s(:,i) = spec_sub_derev_gammatone(s(:,i),fs,H);
% end
%s = WPE_NTT(s);

x{1} = s;
param.fftsize = 2 .^ [10 8]; % frameSize, shift

noise_beg = floor(param.rate*param.tNoiseBeg/param.fftsize(2));
noise_end = floor(param.rate*param.tNoiseEnd/param.fftsize(2));


[X, Xfull, param] = stftAnalyFull(x, param);
[num_mic,~,nFrame,nFreq] = size(Xfull);

t1 = toc;
fprintf('stftAnalyFull duration: %f\n',t1);


PhiNN = zeros(num_mic,num_mic,nFreq);
PhiXX = zeros(num_mic,num_mic,nFreq);


% extract noise covariance matrix
for i = 1: nFreq
    count = 0;
    for j = 1:noise_beg
        PhiNN(:,:,i) = PhiNN(:,:,i) + Xfull(:,:,j,i);
        count = count + 1;
    end
    % for noise period
    for j = nFrame-noise_end+1:nFrame
        PhiNN(:,:,i) = PhiNN(:,:,i) + Xfull(:,:,j,i);
        count = count + 1;
    end
    % normalize
    PhiNN(:,:,i) = PhiNN(:,:,i)./count;
    
    [VN(:,:,i),DN(:,:,i)] = eig(PhiNN(:,:,i));
    % pick the comlumn which has largest eigen value
    VN_vec(:,i) = VN(:,find(diag(DN(:,:,i))==max(diag(DN(:,:,i)))),i);
    count = 0;
    
    % for speech period
    for j = noise_beg + 1:nFrame-noise_end
        PhiXX(:,:,i) = PhiXX(:,:,i) + Xfull(:,:,j,i);
        count = count + 1;
    end
    
    % normalize
    PhiXX(:,:,i) = PhiXX(:,:,i)./count;
    
    if strcmp(subnoisecov,'yes')
        % remove noise covariance matrix from speech covariance matrix
        PhiXX(:,:,i) = PhiXX(:,:,i) - PhiNN(:,:,i);
    end
    
 
    % extract the eigen vector
    [VX(:,:,i),DX(:,:,i)] = eig(PhiXX(:,:,i));
    % pick the comlumn which has largest eigen value
    VX_vec(:,i) = VX(:,find(diag(DX(:,:,i))==max(diag(DX(:,:,i)))),i);
   
end
t2 = toc;
fprintf('Calculate VX_vec, VN_vec duration: %f\n',t2 - t1);
save('PhiNN.mat','PhiNN');
% extract speech noise in the one utterance
% X = num_mic x nFrame x nFreq; VX_vec = num_mic x nFreq; PhiNN = num_mic x num_mic x nFreq


speech_noise_MVDR(1,:,:) = MVDRbeamf(VX_vec,PhiNN,X); 
speech_noise_MVDR(2,:,:) = MVDRbeamf(VN_vec,PhiXX,X);

t3 = toc;
fprintf('MVDR duration: %f\n',t3 - t2);




[A, Y_MVDR] = scaling(X, speech_noise_MVDR, param);

t4 = toc;
fprintf('scaling duration: %f\n',t4 - t3);

yt_mvdr = stftSynth(Y_MVDR(1:2,:,:,:), param);
mvdr_out = yt_mvdr{1}(:,1);

t5 = toc;
fprintf('stftSynth duration: %f\n',t5 - t4);
% audiowrite('speech_mvdr.wav',yt_mvdr{1}(:,1),param.rate);
% audiowrite('noise_mvdr.wav',yt_mvdr{2}(:,1),param.rate);
% if strcmp(derev,'yes')
% end
end

function speech_MVDR = appyling_beamforming(Xfull,PhiNN_mean,X)

[num_mic,~,nFrame,nFreq] = size(Xfull);

for j = 1:nFreq
    PhiXX(:,:,j) = mean(Xfull(:,:,:,j),3);
    PhiXX(:,:,j) = PhiXX(:,:,j) - PhiNN_mean(:,:,j);
    [VX(:,:,j),DX(:,:,j)] = eig(PhiXX(:,:,j));
    VX_vec(:,j) = VX(:,find(diag(DX(:,:,j))==max(diag(DX(:,:,j)))),j);
end
speech_MVDR = MVDRbeamf(VX_vec,PhiNN_mean,X(:,:,:));
end


