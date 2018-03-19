function ret = SaveAllFigures()
% This function allows you to quickly save all currently open figures with
% a custom filename for each in multiple formats.  To use the function
% simply call savefigs with no arguments, then follow the prompts
%
% Upon execution this function will one-by-one bring each currently open
% figure to the foreground.  Then it will supply a text prompt in the main
% console window asking you for a filename.  It will save that figure to
% that filename in the .fig, .emf, .png, and .eps formats.  
%
% The formats that it saves in can be changed by commenting out or adding
% lines below.
%
% Copyright 2010 Matthew Guidry 
% matt.guidry ATT gmail DOTT com  (Email reformatted for anti-spam)

hfigs = get(0, 'children')                          %Get list of figures

hfigsSorted = sortfighandlearray(hfigs,'ascend') 

hfigs = hfigsSorted;

for m = 1:length(hfigs)
    
    % Get figure name
    %filename = get(figure(hfigs(m)),'name');
    %filename(regexp(filename,'[:]')) = [];
    
    filename = 'Plot'
    filename = [filename sprintf('-%s',int2str(m-1))]
    
    
    fpath = 'Figures/';
    
    if strcmp(filename, '0')                        %Skip figure when user types 0
        continue
    else
        saveas(hfigs(m), [fpath filename '.fig']) %Matlab .FIG file
        saveas(hfigs(m), [fpath filename '.eps']) %Windows Enhanced Meta-File (best for powerpoints)
        saveas(hfigs(m), [fpath filename '.png']) %Standard PNG graphics file (best for web)
        
        %eval(['print -depsc2 ' filename])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)
    end
end


end



function [hFigsSorted] = sortfighandlearray(hFigs,varargin)
%SORTFIGHANDLEARRAY Due to changes in Matlab 2014b graphics system, figure
%handles are no longer doubles, but rather an graphics object. Hence, these
%objects can not be sorted directly. This function accepts an array of
%figure handles, and returns an array of similar length, with figure
%handles sorted with respect to their numeric property 'Number'.
%
% *** Inputs:
%   hFigs - Figure handles array. 
%           An array who holds the handles to a series of figures.
%
%   varargin{1} - String.
%           Should match 'ascend' or 'descend', depending on  how the
%           figure handles should be sorted. By default, figure handles are
%           sorted in an ascending order.
%
% *** Outputs:
%   hFigSorted- Figure handles array.
%           Figure handles are sorted with respect to their property
%           'Number'.
% ------------------------------------------------------------------------

narginchk(1,2)

if nargin == 1
   sortStr = 'ascend';
else
   if strcmpi(varargin{1},'ascend') || strcmpi(varargin{1},'descend')
      sortStr = varargin{1};
   else
      error('Bad input. 2nd input argument should match string "ascend" or "descend"')
   end
end

% Allocate output
nFigs = length(hFigs);
hFigsSorted = gobjects(nFigs,1);

% Capture an array of figure numbers
figNumbersArray = [hFigs.Number];

% Sort and find indices of the sorted figures
[~,I] = sort(figNumbersArray,sortStr);

% Assemble output array
hFigsSorted(1:nFigs) = hFigs(I);

end