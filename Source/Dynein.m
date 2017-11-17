classdef Dynein
    properties
        AtpSites = [0, 0] % [primary (hydrolyzable) sites, secondary sites]
        Position = 0;
    end % properties
    properties (Constant)
        NUM_HYDRO_SITES = 1;
        NUM_SITES = 4;
    end
    properties (Dependent)
        S
        PrimaryS
        SecondaryS
    end % dependent properties
    
    methods
        function val = get.S(obj)
            val = sum(obj.AtpSites);
        end % get.S(obj)
        
        % PrimaryS Get/Set
        function val = get.PrimaryS(obj)
            val = obj.AtpSites(1);
        end
        function obj = set.PrimaryS(obj, val)
            obj.AtpSites(1) = val;
        end
        
        % SecondaryS Get/Set
        function val = get.SecondaryS(obj)
            val = obj.AtpSites(2);
        end
        function obj = set.SecondaryS(obj, val)
            obj.AtpSites(2) = val;
        end
        
        function obj = set.AtpSites(obj, val)
            disp(val)
             if length(val) ~= 2
                error('AtpSites should be a vector of length 2.');
            end
            if val(1) > obj.NUM_HYDRO_SITES
                warning(['The number of bound ATPs at primary sites cannot'
                         ' exceed NUM_HYDRO_SITES']);
                val(1) = obj.NUM_HYDRO_SITES;
            end
            if sum(val) > obj.NUM_SITES
                warning(['The total number of bound ATPs cannot exceed'
                         ' NUM_SITES']);
                val(2) = obj.NUM_SITES - obj.NUM_HYDRO_SITES;
            end
            if ~all(val >= 0)
                warning(['The number of bound ATPs cannot be less than'
                         ' 0']);
                val(val < 0) = 0;
            end
            disp(val)
            obj.AtpSites = val;
        end % sel.AtpSites(obj)
        
        function obj = bindAtp(obj)
            if obj.PrimaryS < obj.NUM_HYDRO_SITES
                obj.PrimaryS = obj.PrimaryS + 1;
            elseif obj.S < obj.NUM_SITES
                obj.SecondaryS = obj.SecondaryS + 1;
            else
                warning('Sites full, cannot bind another ATP');
            end
        end
        function obj = unbindAtp(obj)
            if obj.SecondaryS > 0
                obj.SecondaryS = obj.SecondaryS - 1;
            elseif obj.PrimaryS > 0
                obj.PrimaryS = obj.PrimaryS - 1;
            else
                warning('Sites empty, cannot unbind another ATP');
            end
        end
        function obj = hydrolyzeAtp(obj)
            if obj.PrimaryS > 0
                obj.PrimaryS = obj.PrimaryS - 1;
            else
                warning('Primary Sites Empty, cannot hydrolyze ATP');
            end
        end
    end % methods
end % Dynein class