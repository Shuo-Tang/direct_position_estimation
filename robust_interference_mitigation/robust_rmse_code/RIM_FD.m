function x_rim = RIM_FD(x, Th)
%RIM Algorithm
x_rim = zeros(size(x));
for i = 1:size(x, 1)
    
     %% Robustness pre-correlation
     blksize = length(x(i, :));
     fftrawSignal= fft(x(i, :));
%      fftrawSignal= x(i, :);
     MAR = median(abs(fftrawSignal - median(fftrawSignal)));
     sigama=MAR/0.6745; 
%      sigama=std(fftrawSignal);
     sigama = sigama./sqrt(2);
     fftrawSignal_post=zeros(1,blksize);                     
     kF=Th*sigama;%1.345*sigama;
     try
         fftrawSignal_post(abs(fftrawSignal)<=kF)=fftrawSignal(abs(fftrawSignal)<=kF);
         fftrawSignal_post(abs(fftrawSignal)>kF)=kF*fftrawSignal(abs(fftrawSignal)>kF)./abs(fftrawSignal(abs(fftrawSignal)>kF));
    catch
      save("fftrawSignal_post.mat", 'fftrawSignal_post', 'kF', 'fftrawSignal')
      fftrawSignal_post(abs(fftrawSignal)<=kF)=fftrawSignal(abs(fftrawSignal)<=kF);
      fftrawSignal_post(abs(fftrawSignal)>kF)=kF*fftrawSignal(abs(fftrawSignal)>kF)./abs(fftrawSignal(abs(fftrawSignal)>kF));
    end
     
     x_rim(i, :)= ifft(fftrawSignal_post);
%      x_rim(i, :)= fftrawSignal_post;

end
end