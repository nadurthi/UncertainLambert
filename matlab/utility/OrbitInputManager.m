classdef OrbitInputManager < handle
    % [a ecc, inc, W, w, M]
    % Dep_orb_elem(1);%  (a) km
    % Dep_orb_elem(2);%  (ecc) eccentricity
    % Dep_orb_elem(3); % (inc) Inclination (rad)
    % Dep_orb_elem(4);%  (W) Right Ascension of Ascending Node (rad)
    % Dep_orb_elem(5);%  (w) Argument of Perigee (rad)
    % Dep_orb_elem(6);%  (M) Mean Anomaly (rad)
    %
    % % Arrival Orbit Elements
    % Arr_orb_elem(1);% km
    % Arr_orb_elem(2);% eccentricity
    % Arr_orb_elem(3);% Inclination (rad)
    % Arr_orb_elem(4);% Right Ascension of Ascending Node (rad)
    % Arr_orb_elem(5);% Argument of Perigee (rad)
    % Arr_orb_elem(6) ;% Mean Anomaly (rad)
    
    properties
        tle_strings={};
        tle;
        tle_is_set=0;
        a;
        ecc;
        inc;
        W;
        w;
        M;
        orb_elem_vec=[];
        epochDM;  % epoch Date manager
        constants;
    end
    methods
        % constructor. save a copy constants
        function obj= OrbitInputManager(constants)
            obj.constants=constants;
        end
        
        function obj = set_TLE(obj,ss,config)
            obj.tle_strings=ss;
            if length(obj.tle_strings)~=3
                error('TLE needs a cell of length 3');
            end
            obj.tle=tle_cells2orbelem(obj.tle_strings,obj.constants,config);
            obj.epochDM=obj.tle.epochDM;
            obj.tle_is_set=1;
        end
        function obj = set_orb_elements(a,ecc,inc,W,w,M,epochDM)
            obj.a=a;
            obj.ecc=ecc;
            obj.inc=inc;
            obj.W=W;
            obj.w=w;
            obj.epochDM=epochDM;
            obj.M=M;
            
            obj.orb_elem_vec=[a ecc, inc, W, w, M];
        end
        function a=get_a(obj)
            if obj.tle_is_set==1
                a=obj.tle.a;
            else
                a=obj.a;
            end
            
        end
        
        function ecc=get_ecc(obj)
            if obj.tle_is_set==1
                ecc=obj.tle.ecc;
            else
                ecc=obj.ecc;
            end
        end
        
        function inc=get_inc(obj)
            if obj.tle_is_set==1
                inc=obj.tle.inc;
            else
                inc=obj.inc;
            end
        end
        
        function W=get_W(obj)
            if obj.tle_is_set==1
                W=obj.tle.W;
            else
                W=obj.W;
            end
        end
        
        function w=get_w(obj)
            if obj.tle_is_set==1
                w=obj.tle.w;
            else
                w=obj.w;
            end
        end
        function M=get_M(obj)
            if obj.tle_is_set==1
                M=obj.tle.M;
            else
                M=obj.M;
            end
        end
        function T=get_period(obj)
            if obj.tle_is_set==1
                T=obj.tle.T;
            else
                T=2*pi/sqrt(obj.constants.mu/obj.get_a()^3);
            end
        end
        function orb_elem_vec=get_orb_elem_vec(obj)
            % [a ecc, inc, W, w, M]
            orb_elem_vec=[obj.get_a(), obj.get_ecc(),obj.get_inc(),obj.get_W(),obj.get_w(),obj.get_M()];
        end
        function [r,v]=get_rv(obj)
           [r,v]=elm2rv(obj.get_a(), obj.get_ecc(),obj.get_inc(),obj.get_W(),obj.get_w(),obj.get_M(),0,obj.constants.mu); 
        end
        function ss=get_as_cell_string(obj)
           ss={['a : ',num2str(obj.get_a())],...
               ['e : ',num2str(obj.get_ecc()) ],...
               ['i : ',num2str(obj.get_inc()) ],...
               ['w : ',num2str(obj.get_w()) ],...
               ['W : ',num2str(obj.get_W()) ],...
               ['M : ',num2str(obj.get_M()) ],...
               ['T : ',num2str(obj.get_period()) ],...
               ['epoch : ',obj.epochDM.date_string]} ;
        end
        function print_summary(obj)
           ss=obj.get_as_cell_string();
           disp(ss');
        end
    end
end