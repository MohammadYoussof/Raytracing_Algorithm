function [h, t, gather] = SEGY_ReadCMPGather(segy, cmpnum, binsize)

% [h, t, gather] = SEGY_ReadCMPGather(segy, cmpnum, binsize)
%
% Returns the CMP gather binned with binsize at segy.cmps(cmpnum)
%
% Chad Hogan, 2008
%
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

% $Id: SEGY_ReadCMPGather.m,v 1.1 2008/03/04 22:37:12 cmhogan Exp $

if isnan(segy.cmps)
    error('You did not find CMPs! Use SEGY_FindCMPs() to find them.');
end

if (cmpnum > length(segy.cmps))
    error('Aint that many CMPs gathers, hotshot.');
end
bin = binsize/2;

disp(['Going to search ' num2str(segy.numtraces) ' traces']);

gsize = 0;

cmp = segy.cmps(cmpnum);

for idx = 1:segy.numtraces
    SEGY_TraceSeek(segy, idx);
    fseek(segy.FILE, 72, 0);    % move to sx
    thissx = fread(segy.FILE, 1, 'int');
    thissy = fread(segy.FILE, 1, 'int');
    thisgx = fread(segy.FILE, 1, 'int');
    thisgy = fread(segy.FILE, 1, 'int');

    cmpx = (thissx + thisgx) / 2;
    cmpy = (thissy + thisgy) / 2;
    
    if (idx == 1)
        x0 = cmpx;
        y0 = cmpy;
    end
    
    dx = cmpx - x0;
    dy = cmpy - y0;
    
    dist = sqrt(dx^2 + dy^2);
    
    dist = round(dist / bin) * bin;
    
    if(dist == cmp)
        trace = SEGY_ReadTrace(segy, idx);
        gather(:, gsize + 1) = trace.data;
        thish = sqrt((thisgx - thissx).^2 + (thisgy - thissy).^2);
        if(thisgx == thissx)
            thish = sign(thisgy - thissy) * thish;
        else
            thish = sign(thisgx - thissx) * thish;
        end
        h(gsize+1) = thish;
        gsize = gsize + 1;
    end
    
    if(mod(idx, 1000) == 0)
        disp(['done ' num2str(idx) ' of ' num2str(segy.numtraces)]);
    end
end

t = ((1:size(gather, 1))-1) *segy.bhead.hdt;