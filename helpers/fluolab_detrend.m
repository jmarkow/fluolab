function [DATA,BASELINE]=fluolab_detrend(DATA,varargin)
%simple detrending routine
%
%
%
%
%

% parameter collection

nparams=length(varargin);
dff=1; % convert to dff or just detrend?
fs=22; % sampling rate of signal
win=.4; % size of window (to the left and right of each point, hence .4 is a .8 window in total length)
per=8; % used for the prctile method, which prctile to use
method='prctile'; % 'prctile','regression'
sliding=1;
edge='r';

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'per'
			per=varargin{i+1};
		case 'dff'
			dff=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'win'
			win=varargin{i+1};
		case 'sliding'
			sliding=varargin{i+1};
		case 'edge'
			edge=varargin{i+1};
	end
end

% ensure proper formatting

if isvector(DATA), DATA=DATA(:); end

win_samples=round(win*fs);

[nsamples,ntrials]=size(DATA);

win_len=length([-win_samples:win_samples]);
win_x=[1:win_len]';
win_intercept=ones(size(win_x));

NEWDATA=DATA;
BASELINE=zeros(size(NEWDATA));

if win>0
	for i=1:ntrials

		curr_data=DATA(:,i);

		% repeat data at the edges

		curr_data=[ repmat(curr_data(1),[win_samples 1]);curr_data;repmat(curr_data(end),[win_samples 1]) ];

		counter=1;

		% slide the window

		for j=win_samples+1:nsamples+win_samples

			idx=j-win_samples:j+win_samples;
			tmp=curr_data(idx);

			switch lower(method(1))

				case 'p'

					% sliding percentile

					tmp_baseline=prctile(tmp,per);

					if dff
						tmp=((DATA(j-win_samples,i)-tmp_baseline)./tmp_baseline)*100;
					else
						tmp=(DATA(j-win_samples,i)-tmp_baseline);
					end

				case 'r'

					% least squares solution, brute force

					tmp=tmp(:);

					tmp_coeffs=[win_intercept win_x]\tmp;
					tmp_baseline=win_x.*tmp_coeffs(2)+tmp_coeffs(1);

					if dff
						tmp=((tmp-tmp_baseline)./tmp_baseline)*100;
					else
						tmp=(tmp-tmp_baseline);
					end

					tmp=tmp(win_samples+1);

				otherwise

					error('Did not recognize detrend method...')

			end

				% replace sample with detrended version
			NEWDATA(j-win_samples,i)=tmp;
			BASELINE(j-win_samples,i)=tmp_baseline;

		end

	end
else

	switch lower(method(1))
		case 'p'

			baseline=repmat(prctile(NEWDATA,per),[nsamples 1]); % baseline per trial

			if dff
				NEWDATA=((NEWDATA-baseline)./baseline)*100;
			else
				NEWDATA=NEWDATA-baseline;
			end

			BASELINE=baseline;

		case 'r'

			win_x=[1:nsamples]';
			win_intercept=ones(size(win_x));

			for i=1:ntrials
				tmp=NEWDATA(:,i);
				tmp=tmp(:);

				tmp_coeffs=[win_intercept win_x]\tmp;
				tmp_baseline=win_x*tmp_coeffs(2)+tmp_coeffs(1);

				if dff
					tmp=((tmp-tmp_baseline)./tmp_baseline)*100;
				else
					tmp=(tmp-tmp_baseline);
				end

				NEWDATA(:,i)=tmp;
				BASELINE(:,i)=tmp_baseline;
			end
	end
end

DATA=NEWDATA; % overwrite old data with detrended data
clear NEWDATA;
