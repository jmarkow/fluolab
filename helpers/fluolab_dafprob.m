function [TIMEDIFF,COUNT,PROB,MIN_T]=fluolab_dafprob(DATENUMS,TRIALS,varargin)
%
%
%
%
%

history=60*10;


% get all to all time difference in seconds


datevecs=datevec(DATENUMS);
ndates=length(DATENUMS);
TIMEDIFF=zeros(ndates,ndates);

for i=1:ndates
	for j=1:ndates
		TIMEDIFF(i,j)=etime(datevecs(i,:),datevecs(j,:));
	end
end

% for each timebin, how many daf trials in past 30 minutes
%

PROB=zeros(1,ndates);
MIN_T=zeros(1,ndates);
COUNT=zeros(1,ndates);

for i=1:ndates
	idx=find(TIMEDIFF(i,:)>0&TIMEDIFF(i,:)<history);
	ndaf=length(intersect(TRIALS.all.daf,idx));
	PROB(i)=ndaf/(length(idx)+eps);
    COUNT(i)=ndaf;
    tmp=min(TIMEDIFF(i,TIMEDIFF(i,:)>0));
    if ~isempty(tmp)
        MIN_T(i)=tmp;
    end
end




