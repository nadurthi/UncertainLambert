classdef DateManager  < handle
    % DataManager handles all the conversions between formats
    % SPICE ET(Ephemeris Time) is included.
    % Optionally, one can use Matlab inbuilt datenum to convert dates
    % formats
    
    % Use as
    % '2017-11-19 14:16:57'
    % DM=DateManager()
    % DM.set_date_string('2017-11-19 14:16:57')
    % DM.get_ET() % get SPICE Ephemeris Time (ET)
    
    properties
        config;
        
        date_str_format='yyyy-mm-dd HH:MM:SS';
        spice_date_str_format='YYYY-MM-DD HR:MN:SC';
        date_string
        date_ET=NaN; %Ephemeris Time
        
        date_Y
        date_dayfrac
        date_indays_from_Y01M00D00
        date_inhours_from_Y01M00D00
        date_insecs_from_Y01M00D00
    end
    methods
        function obj =DateManager(config)
            obj.config=config;
        end
        function obj = set_date_year_dayfrac(obj,Y,dayfrac)
            obj.date_Y=Y;
            obj.date_dayfrac=dayfrac;
            
            obj.date_indays_from_Y01M00D00=datenum(Y,0,dayfrac);
            obj.date_inhours_from_Y01M00D00=obj.date_indays_from_Y01M00D00*24;
            obj.date_insecs_from_Y01M00D00=obj.date_inhours_from_Y01M00D00*3600;
            obj.date_string=datestr( obj.date_indays_from_Y01M00D00 ,obj.date_str_format);
        end
        
        function obj= set_date_string(obj,ss)
            dtvec=datevec(ss,obj.date_str_format);
            obj.date_string=ss;
            obj.date_indays_from_Y01M00D00=datenum(ss,obj.date_str_format);
            obj.date_inhours_from_Y01M00D00=obj.date_indays_from_Y01M00D00*24;
            obj.date_insecs_from_Y01M00D00=obj.date_inhours_from_Y01M00D00*3600;
            
            obj.date_Y=dtvec(1);
            obj.date_dayfrac=obj.date_indays_from_Y01M00D00-datenum(obj.date_Y,0,0);
            
        end
        function obj=set_date_insecs_from_Y01M00D00(obj,insecs)
            
            obj.date_insecs_from_Y01M00D00=insecs;
            obj.date_inhours_from_Y01M00D00=obj.date_insecs_from_Y01M00D00/3600;
            obj.date_indays_from_Y01M00D00=obj.date_inhours_from_Y01M00D00/24;
            obj.date_string=datestr(obj.date_indays_from_Y01M00D00,obj.date_str_format);
            dtvec=datevec(obj.date_string,obj.date_str_format);
            obj.date_Y=dtvec(1);
            obj.date_dayfrac=obj.date_indays_from_Y01M00D00-datenum(obj.date_Y,0,0);
            
        end
        
        function obj= set_date_ET(obj,ET)
            if obj.config.useSPICE==1
                ss=cspice_timout(ET,obj.spice_date_str_format);
                obj=obj.set_date_string(ss);
            else
                obj=obj.set_date_insecs_from_Y01M00D00(ET);
            end
        end
        
        function ET = get_ET(obj)
            if obj.config.useSPICE==1
                if isnan(obj.date_ET)
                    string = obj.date_string;
                    obj.date_ET = cspice_str2et( string );
                end
                ET=obj.date_ET;
            else
                ET=obj.date_insecs_from_Y01M00D00;
            end
        end
        
        function ss=get_ET_str(obj,format)
            if obj.config.useSPICE==1
                ss=cspice_timout(obj.date_ET,format);
            else
                
                ss=obj.date_string;
            end
            
        end
    end
end