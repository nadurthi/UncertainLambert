function OE=tle_cells2orbelem(ss,constants,config)

% NOAA 14                 
% 1 23455U 94089A   97320.90946019  .00000140  00000-0  10191-3 0  2621
% 2 23455  99.0090 272.6745 0008546 223.1686 136.8816 14.11711747148495

% GPS BIIR-2  (PRN 13)    
% 1 24876U 97035A   17281.24622173  .00000020  00000-0  00000+0 0  9998
% 2 24876  55.5291 222.6639 0032993  99.2740 261.0932  2.00565211148090

% give input as a matlab cell of length 3:
% ss={  'GPS BIIR-3  (PRN 11)',
%		'1 25933U 99055A   17281.42067461 -.00000010  00000-0  00000+0 0  9990',
%		'2 25933  51.7133  72.6825 0168212  97.6597 316.4246  2.00555860131919'; }



% Two-line element set

mu = constants.mu; %  Standard gravitational parameter for the earth

% 1957 is the year satellites were launched
sat_num=ss{2}(3:7);
classification=ss{2}(8);
int_desg_launchyear=ss{2}(10:11);
if str2double(int_desg_launchyear)>=57
    int_desg_launchyear=str2double(['19',int_desg_launchyear]);
else
   int_desg_launchyear=str2double(['20',int_desg_launchyear]); 
end

int_desg_launchnum=ss{2}(12:14);
int_desg_piecelaunch=ss{2}(15:17);
epoch_year=ss{2}(19:20);
if str2double(epoch_year)>=57
    epoch_year=str2double(['19',epoch_year]);
else
   epoch_year=str2double(['20',epoch_year]); 
end

epoch_dayfrac=str2double(ss{2}(21:32));
% epoch time with respect to '1-jan-1957'
% we might need a much better to deal with the epoch time
% epoch_t=convertepoch2time(epoch_year,epoch_dayfrac);

% remember datenum is number (+fraction) of the number of days from Jan-00-0000 
epochDM=DateManager(config);
epochDM=epochDM.set_date_year_dayfrac(epoch_year,epoch_dayfrac);

% epoch.year=epoch_year;
% epoch.dayfrac_of_year=epoch_dayfrac;

% epoch.datenum_indays=datenum(epoch_year,0,epoch_dayfrac);
% epoch.datenum_inhours=epoch.datenum_indays*24;
% epoch.datenum_insecs=epoch.datenum_inhours*3600;
% epoch.datestr=datestr( epoch.datenum_indays ,'mm-dd-yyyy HH:MM:SS');



bstar_drag=ss{2}(54:61);
if bstar_drag(1)=='-'
    bstar_drag=-str2double(bstar_drag(2:6))*eval(['1e-',bstar_drag(8)])/1e5;
else
    bstar_drag=str2double(bstar_drag(2:6))*eval(['1e-',bstar_drag(8)])/1e5;
end

inc=str2double(ss{3}(9:16))*pi/180;
W=str2double(ss{3}(18:25))*pi/180;
ecc=str2double(ss{3}(27:33))/1e7;
w=str2double(ss{3}(35:42))*pi/180;
M=str2double(ss{3}(44:51))*pi/180;
n=str2double(ss{3}(53:63))*2*pi/(24*3600);

a = (mu/n^2)^(1/3);     % Semi-major axis [km]  

T=2*pi/sqrt(constants.mu/a^3);


(mu/n^2)

% Six orbital elements 
OE.a=a;
OE.ecc=ecc;
OE.inc=inc;
OE.W=W;
OE.w=w;
OE.M=M;
OE.n=n
OE.T=T;
OE.epochDM=epochDM;
% OE.epoch_dayfrac=epoch_dayfrac;

% OE.epoch_datenum=epoch_datenum;
% OE.epoch_datestr=epoch_datestr;

OE.bstar_drag=bstar_drag;
OE.tle=ss;
OE.sat_num=sat_num;
OE.classification=classification;
OE.int_desg_launchyear=int_desg_launchyear;
OE.int_desg_launchnum=int_desg_launchnum;
OE.int_desg_piecelaunch=int_desg_piecelaunch;

% OE = [a e inc W w M E n epoch Db];
% fprintf('\n a [km]   e      inc [deg]  RAAN [deg]  w[deg]    E [deg] \n ')
% fprintf('%4.2f  %4.4f   %4.4f       %4.4f     %4.4f    %4.4f', OE);
