function headerobj=SEGY_StandardizeHeader(headerobj)
% tracehhead=SEGY_StandardizeHeader(tracehead);
% binaryhead=SEGY_StandardizeHeader(binaryhead);
% trace=SEGY_StandardizeHeader(trace);
%
% SEGY_StandardizeHeader is a function that allows the user to edit or
%   create a TextHeader Object.
%
% Inputs:
%   headerobj can either be a BinaryHeader Object, a TraceHeader Object or
%     a Trace Object
%
% Outputs:
%   headerobj will return the same type of Object as given in the inputs.
%
%
% Heather Lloyd 2010, Kevin Hall 2009, Chad Hogan 2004
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
%
if isa(headerobj,'Trace')
    hdrobj=headerobj.traceheader;
    headerobj.traceheader=uiHeaderConverter(hdrobj);
else
    headerobj=uiHeaderConverter(headerobj);
end