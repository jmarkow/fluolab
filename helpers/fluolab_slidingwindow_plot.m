function [STEPS,t]=fluoflab_sliding_window(T,F,S,DATA,DATA_T,WIN,WIN_OVERLAP,varargin)
%
%
%
%
colors='jet';
cbar_dist=.1;
cbar_width = .025; 
datenums=[];
time_dist=.005;
time_width=.075;
escape_dist=.005;
escape_width=.075;
escape_idx=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'cbar_dist'
			cbar_dist=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'cbar_width'
			cbar_width=varargin{i+1};
		case 'datenums'
			datenums=varargin{i+1};
		case 'label'
			label=varargin{i+1};
		case 'escape_idx'
			escape_idx=varargin{i+1};
	end
end

[nsamples,ntrials]=size(DATA);


isdate=0;
isescape=0;
if length(datenums)==ntrials
	isdate=1;
end

if length(escape_idx)==ntrials
	isescape=1;
end

if isescape & isdate
	%time_width=time_width+escape_width;
end

step_size=WIN-WIN_OVERLAP;
STEPS=1:step_size:ntrials-WIN-1;
sliding_fluo=zeros(nsamples,length(STEPS));

counter=1;
for i=STEPS
	sliding_fluo(:,counter)=mean(DATA(:,[i:i+WIN])');
	counter=counter+1;
end


if length(datenums)==ntrials
	isdate=1;
end

if length(escape_idx)==ntrials
	isescape=1;
end


ax=[];

nplots=4;

idx=1;

ax(end+1)=subplot(nplots,1,idx);
imagesc(T,F,S);axis xy;box off;
set(gca,'xtick',[]);
title(['Win ' num2str(WIN) ' overlap ' num2str(WIN_OVERLAP)]);
colormap(colors);

idx=idx+1;

ax(end+1)=subplot(nplots,1,idx:idx+2);
imagesc(DATA_T,[],sliding_fluo');
ylabel('Trials');
xlabel('Time(s)');
box off;
linkaxes(ax,'x');
xlim([max(DATA_T(1),T(1)) min(DATA_T(end),T(end))]);
title([label]);


pos=get(gca,'position');
clims=caxis();
clims=round(clims*100)/100;
caxis(clims);
set(ax(end),'TickLength',[0 0]);
set(ax(end),'Position',[pos(1) pos(2)+cbar_dist+cbar_width pos(3) pos(4)-cbar_dist-cbar_width]);

hc = colorbar('location','southoutside','position',...
	[pos(1) pos(2) pos(3) cbar_width],'fontsize',10,'xtick',clims);

newpos=get(ax(end),'position');


idx=idx+3;
dateflag=0;


if isdate

	counter=1;
	tmp_datenums=zeros(1,length(STEPS));
	for i=STEPS
		tmp_datenums(counter)=mean(datenums(i:i+WIN));
		counter=counter+1;
	end

	if (length(tmp_datenums)<2) | (tmp_datenums(end)<=tmp_datenums(1))
		return;
	end
	
	% rescale all axes prior to insertion

	% show how time progresses across trials

	% make a custom axis with time progression on the right

	for i=1:length(ax)
		pos=get(ax(i),'position');
		set(ax(i),'position',[pos(1) pos(2) pos(3)-time_width-time_dist pos(4)]);
	end

	newpos=get(ax(end),'position');
	pos=get(hc,'position');
	set(hc,'position',[pos(1) pos(2) newpos(3) pos(4)]);

	ax(end+1)=axes('position',[newpos(1)+newpos(3)+time_dist newpos(2) time_width newpos(4)]);

	plot(tmp_datenums,[1:length(tmp_datenums)],'linewidth',1.5)	
	set(gca,'XTick',[tmp_datenums(1) tmp_datenums(end)],'FontSize',11);
	set(gca,'YDir','rev','YTick',[]);
	%set(gca,'YColor',get(gcf,'color'));
	xlim([tmp_datenums(1) tmp_datenums(end)])
	ylim([1 length(tmp_datenums)])

	xtick_labels{1}=datestr(tmp_datenums(1),'HH:MM (PM)');
	xtick_labels{2}=datestr(tmp_datenums(end),'HH:MM (PM)');

	xtick_pos=get(ax(end),'XTick');
	set(ax(end),'xtick',[]);

	% set y shift to 10% axis
	
	axis_range=length(tmp_datenums);

	linkaxes(ax(end-1:end),'y');
	text(xtick_pos,repmat(-.05*axis_range,[length(xtick_pos) 1]),xtick_labels,'rotation',45);
	dateflag=1;

end


if isescape

	counter=1;
	%tmp_escape=zeros(1,length(STEPS));
	tmp_escape=escape_idx(STEPS);

	if (length(tmp_escape)<2) 
		return;
	end

	if isdate
		n=length(ax)-1;
	else
		n=length(ax);
	end

	for i=1:n
		pos=get(ax(i),'position');
		set(ax(i),'position',[pos(1) pos(2) pos(3)-escape_width-escape_dist pos(4)]);
	end

	if isdate
		pos=get(ax(end),'position');
		set(ax(end),'position',[pos(1)-escape_width-escape_dist pos(2) pos(3) pos(4)]);
	end

	newpos=get(ax(1),'position');
	pos=get(hc,'position');
	set(hc,'position',[pos(1) pos(2) newpos(3) pos(4)]);
	newpos=get(ax(end),'position');

	ax(end+1)=axes('position',[newpos(1)+newpos(3)+escape_dist newpos(2) escape_width newpos(4)]);

	plot(tmp_escape,[1:length(tmp_escape)],'linewidth',1.5);
	set(gca,'XTick',[0 1],'FontSize',11);
	set(gca,'YDir','rev','YTick',[]);

	xlim([0 1]);
	ylim([1 length(tmp_escape)]);

	xtick_pos=get(ax(end),'XTick');
	set(ax(end),'xtick',[]);

	axis_range=length(tmp_escape);

	linkaxes(ax(end-1:end),'y');
	%text(xtick_pos,repmat(-.05*axis_range,[length(xtick_pos) 1]),xtick_labels,'rotation',45);

	escapeflag=1;

end



