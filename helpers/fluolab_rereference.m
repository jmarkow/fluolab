function [SIGNAL,BASELINE]=fluolab_rereference(SIGNAL,REFERENCE)
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

opts=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','TolX',1e-20,'TolFun',1e-20,'Display','Off');

coeffs=zeros(2,trials1);

for i=1:trials1
  fun=@(x) mean((((x(1)*REFERENCE(:,i)))-SIGNAL(:,i)).^2);
  coeffs(:,i)=lsqnonlin(fun,[.8 0],[],[],opts);
end

m=repmat(coeffs(1,:),[samples1 1]);
b=repmat(coeffs(2,:),[samples1 1]);
BASELINE=REFERENCE.*m+b;
SIGNAL=(SIGNAL-BASELINE);
