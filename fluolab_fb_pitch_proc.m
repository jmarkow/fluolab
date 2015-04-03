function [pitch]=fluolab_pitch_proc(AUDIO,TARGET,varargin)
%
%
%
%
%
%
%

if nargin<2, error('Insufficient number of arguments'); end

cf=[1e3:1e3:6e3];
order=201;
stopband=repmat(250,[1 length(cf)]);
bw=repmat(150,[1 length(cf)]);

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'order'
			order=varargin{i+1};
		case 'cf'
			cf=varargin{i+1};
		case 'stopband'
			stopband=varargin{i+1};
		case 'bw'
			bw=varargin{i+1};
	end
end

%% build filterbank

[filterbank,filterbank_resp]=sylldet_filterbank(AUDIO.fs,'cf',cf,'stopband',stopband,'bw',bw,'order',order);

%% test responses

[nsamples,ntrials]=size(AUDIO);

nfilters=length(filterbank);

pitch.val=[];
pitch.t=[];

for i=1:nfilters
	response=filter(filterbank(i),double(AUDIO.data.^2));
	[pitch.val(:,:,i),pitch.t]=sylldet_pitch_batch(response,AUDIO.fs,'len',5,'overlap',0,'filtering',[]);
	pitch.val(:,:,i)=pitch.val(:,:,i)/i;
end

% collect pitch values from the target

[~,idx1]=min(abs(pitch.t-TARGET(1)))
[~,idx2]=min(abs(pitch.t-TARGET(2)))

%

pitch.target.mat=pitch.val(idx1:idx2,:,:);
pitch.target.mu.trial=mean(pitch.target.mat,3); % average over samples and filters
pitch.target.mu.trial_sample=mean(pitch.target.mu.trial);

% collect pitch values
