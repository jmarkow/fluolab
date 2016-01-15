function SIGNAL=fluolab_rereference(SIGNAL,REFERENCE)
%
%
%
%
%

% check dimensions

[samples1,trials1]=size(SIGNAL);
[samples2,trials2]=size(REFERENCE);

if samples1~=samples2
  error('Unequal samples for rereference');
end

if trials1~=trials2
  error('Unequal trial number for rereference');
end

opts=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display','off');

coeffs=zeros(2,trials1);

for i=1:trials1
  fun=@(x) mean((((x(1)*REFERENCE(:,i))+x(2))-SIGNAL(:,i)).^2);
  coeffs(:,i)=lsqnonlin(fun,[.8 .8],[],[],opts);
end

m=repmat(coeffs(1,:),[samples1 1]);
b=repmat(coeffs(2,:),[samples1 1]);

baseline=REFERENCE.*m+b;
SIGNAL=(SIGNAL-baseline)./baseline;
