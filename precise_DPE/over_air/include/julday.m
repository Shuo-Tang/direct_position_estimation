function jd = julday(y,m,d,h)

% SYNTAX:
%   jd = julday(y,m,d,h);
%
% INPUT:
%   y = year
%   m = month
%   d = day
%   h = hour
%
% OUTPUT:
%   jd = julian day
%
% DESCRIPTION:
%   Julian day computation.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) Kai Borre
% Kai Borre 02-14-01
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

if m <= 2
    y = y-1;
    m = m+12;
end

%return julian day
jd = floor(365.25*(y+4716))+floor(30.6001*(m+1))+d+h/24-1537.5;


%% GPSTk converter

%    long convertCalendarToJD( int yy, 
%                              int mm,
%                              int dd ) 
%       throw()
%    {
%       if(yy == 0)
%          --yy;         // there is no year 0
% 
%       if(yy < 0) 
%          ++yy;
%       
%       long jd;
%       double y = static_cast<double>( yy ), 
%          m = static_cast<double>( mm ), 
%          d = static_cast<double>( dd );
% 
%          // In the conversion from the Julian Calendar to the Gregorian
%          // Calendar the day after October 4, 1582 was October 15, 1582.
%          //
%          // if the date is before October 15, 1582
%       if(yy < 1582 || (yy == 1582 && (mm < 10 || (mm == 10 && dd < 15))))
%       {
%          jd = 1729777 + dd + 367 * yy 
%             - static_cast<long>(7 * ( y + 5001 +
%                                       static_cast<long>((m - 9) / 7)) / 4) 
%             + static_cast<long>(275 * m / 9);
%       }
%       else   // after Oct 4, 1582
%       {     
%         jd = 1721029 + dd + 367 * yy 
%            - static_cast<long>(7 * (y + static_cast<long>((m + 9) / 12)) / 4)
%            - static_cast<long>(3 * (static_cast<long>((y + (m - 9) / 7) / 100) 
%                                     + 1) / 4) 
%            + static_cast<long>(275 * m / 9);
% 
%             // catch century/non-400 non-leap years
%          if( (! (yy % 100) && 
%               (yy % 400) && 
%               mm > 2 && 
%               mm < 9)      || 
%              (!((yy - 1) % 100) &&
%               ((yy - 1) % 400) &&
%               mm == 1)) 
%          {
%             --jd;
%          }
%       }
%       return jd;
%    }
%% Time to SOD
%    double convertTimeToSOD( int hh, 
%                             int mm,
%                             double sec ) 
%       throw()
%    {
%       return (sec + 60. * (mm + 60. * hh));
%    }