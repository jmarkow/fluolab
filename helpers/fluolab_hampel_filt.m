function [TO_DEL]=fluolab_hampel_filt(DATA,varargin)
% check how flat the spectrum is each timepoint 
% regression?

% the template cutoff could be defined by the 95th prctile of the abs(noise) magnitude

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

nmads=6;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'nmads'
			nmads=varargin{i+1};
	end
end

[nsamples,ntrials]=size(DATA);

med=repmat(median(DATA,2),[1 ntrials]);
diffmad=DATA-med;

posmad=zeros(nsamples,1);
negmad=zeros(nsamples,1);

for i=1:nsamples
	posmad(i)=abs(median(diffmad(i,diffmad(i,:)>0)));
	negmad(i)=abs(median(diffmad(i,diffmad(i,:)<0)));
end

%posmad=abs(median(diffmad(diffmad>0)))
%negmad=abs(median(diffmad(diffmad<0)))

posmad=repmat(posmad,[1 ntrials]);
negmad=repmat(negmad,[1 ntrials]);

[~,TO_DEL]=find((DATA>(med+nmads*posmad))|(DATA<(med-nmads*negmad)));

% take the song/nonsong power ratio
