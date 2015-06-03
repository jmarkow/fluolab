function [DATA,MIC_DATA,DIST]=ephys_collectsilence(DIR,varargin)
%
%
%
%

% DIR specifies directory to process, GAP specifies how long 

if nargin<1
	error('ephysPipeline:tfhistogram:notenoughparams','Need 1 argument to continue, see documentation');
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

min_dist=1; % minimum distance from any vocalization (in seconds,typically 4 seconds)
seg_length=1; % how long should the segments be (typically 1 second)
max_trials=inf; % maximum number of trials to collect

song_ratio=2; % power ratio between song and non-song band
song_len=.005; % window to calculate ratio in (ms)
song_overlap=0; % just do no overlap, faster
song_thresh=.3; % between .2 and .3 seems to work best (higher is more exclusive)
song_band=[2e3 6e3];
song_pow=-inf; % raw power threshold (so extremely weak signals are excluded)
song_duration=.8; % moving average of ratio
filtering=300; % changed to 100 from 700 as a more sensible default, leave empty to filter later
audio_pad=.05;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'min_dist'
			min_dist=varargin{i+1};
		case 'seg_length'
			seg_length=varargin{i+1};
		case 'savedir'
			savedir=varargin{i+1};
		case 'max_trials'
			max_trials=varargin{i+1};
		case 'song_thresh'
			song_thresh=varargin{i+1};
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) cycle through each file in the directory
% 2) for each file run song detection
% 3) collect all seg_length long samples that are min_dist from any vocalizations
% 4) aggregate and save for analysis

% use typical song detection, power in the relevant frequency bands
% first get the file list

filelist=dir(fullfile(DIR,'*.mat'));

DATA=[];
MIC_DATA=[];
DIST=[];

counter=1;
contflag=1;

length(filelist)
for i=1:length(filelist)

	i
	load(fullfile(DIR,filelist(i).name),'adc','audio');

	[b,a]=butter(5,[filtering/(audio.fs/2)],'high');
	audio.norm_data=filtfilt(b,a,double(audio.data));
	audio.norm_data=audio.norm_data./max(abs(audio.norm_data));

	% filter the data then run through song_det with extremely low threshold	

	[song_bin,song_t]=zftftb_song_det(audio.norm_data,audio.fs,'song_band',song_band,...
				'len',song_len,'overlap',song_overlap,'song_duration',song_duration,...
				'ratio_thresh',song_ratio,'song_thresh',song_thresh,'pow_thresh',song_pow);

	raw_t=[1:length(audio.norm_data)]./audio.fs;

	% interpolate song detection to original space, collate idxs

	detection=interp1(song_t,double(song_bin),raw_t,'nearest'); 
	detection=detection>0;
	detection(isnan(detection))=0;
	detection_inv=~detection;
	
	ext_pts_silence=markolab_collate_idxs(detection_inv,round(audio_pad*audio.fs))
	ext_pts_song=markolab_collate_idxs(detection,round(audio_pad*audio.fs));

	% where was song detected (song_bin>0)

	for j=1:size(ext_pts_silence,1)

		% trailing edge of the song epoch
		
		len_smps=diff(ext_pts_silence(j,:));
		len_t=len_smps/audio.fs;
		seg_smps=round(seg_length*audio.fs);
		dist_smps=round(min_dist*audio.fs);

		start_point=ext_pts_silence(j,1);
		stop_point=ext_pts_silence(j,2);

		% sliding window starting at start_point+dist_smps

		left_edge=start_point;
		step_size=seg_smps;
		right_edge=left_edge+step_size;

		while right_edge<stop_point

			song_dist_left=min(abs(left_edge-ext_pts_song(:)))
			song_dist_right=min(abs(right_edge-ext_pts_song(:)))

			if isempty(ext_pts_song)
				song_dist_left=inf;
				song_dist_right=inf;
			end

			flag_left=song_dist_left>dist_smps;
			flag_right=song_dist_right>dist_smps;
			
			% check candidate
			
			[b,a]=ellip(3,.2,40,[3e3 7e3]/(audio.fs/2),'bandpass');

			mic_amp=filtfilt(b,a,audio.data(left_edge:right_edge));
			mic_amp=20*log10(sqrt(markolab_smooth(mic_amp.^2,round(.01*24.414e3))));
			flag_amp=mean(mic_amp)<-60;

			if flag_left & flag_right & flag_amp
				DATA(:,counter)=adc.data(left_edge:right_edge,1);
				MIC_DATA(:,counter)=audio.data(left_edge:right_edge);
				DIST(:,counter)=[song_dist_left/audio.fs;song_dist_right/audio.fs];
				counter=counter+1
			end
			
			left_edge=left_edge+seg_smps;
			right_edge=right_edge+seg_smps;

		end	

		if counter>max_trials
			contflag=0;
		end

	end

	if ~contflag
		break;
	end

end

% delete unused 

