function change_tickLabels(varargin)

textopts = {};
if length(varargin{1})==1 && ...
        ishandle(varargin{1}) && ...
        strcmpi(get(varargin{1},'Type'),'axes');
    Ha = varargin{1};
    ytickpos = varargin{2};
    ytickstring = varargin{3};
    if nargin > 3
        textopts = varargin(4:end);
    end
else
    Ha = gca;
    Hfig = get(Ha,'Parent');
    ytickpos = varargin{1};
    ytickstring = varargin{2};
    if nargin > 2
        textopts = varargin(3:end);
    end
end

set(Ha,'Ytick',ytickpos, 'Yticklabel','')
h_olds = findobj(Ha, 'Tag', 'MUXTL');
if ~isempty(h_olds)
    delete(h_olds)
end

%% Make Yticklabels
NTick = length(ytickpos);
Ybot = min(get(gca,'YLim'));
ht = zeros(NTick,1);
for ii = 1:NTick
    ht(ii) = text('String',ytickstring{ii}, ...
        'Units','data', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center ', ...
        'Position',[ytickpos(ii) Ybot], ...
        'Tag','MUXTL');
end
if ~isempty(textopts)
    set(ht,textopts{:})
end

%% squeeze axis if needed

set(Ha,'Units','pixels')
Axpos = get(Ha,'Position');
% set(Hfig,'Units','pixels')
% Figpos = get(Hfig,'Position');

set(ht,'Units','pixels')
TickExt = zeros(NTick,4);
for ii = 1:NTick
    TickExt(ii,:) = get(ht(ii),'Extent');
end

needmove = -(Axpos(2) + min(TickExt(:,2)));

if needmove>0;
    Axpos(2) = Axpos(2)+needmove+2;
    Axpos(4) = Axpos(4)-needmove+2;
    set(Ha,'Position',Axpos);
end

set(Ha,'Units','normalized')
set(ht,'Units','normalized') 