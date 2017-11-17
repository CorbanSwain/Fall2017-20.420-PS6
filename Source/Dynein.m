classdef Dynein
    properties
        Position = 0;
    end % properties
    properties (SetAccess = 'private')
        AtpSites = [0, 0] % [primary (hydrolyzable) sites, secondary sites]
    end
    properties (Constant)
        NUM_HYDRO_SITES = 1
        NUM_SITES = 4
        STEP_INCREMENT = 8
        K_ON_1_0 = 4E5 % 1/M 1/s
        K_OFF_1_0 = 10 % 1/s
        K_OFF_2_0 = 250 % 1/s
        K_CAT = 55 % 1/s
        K_ON_0 = [1, 1, 1/4, 1/6] * Dynein.K_ON_1_0
        K_OFF_0 = [Dynein.K_OFF_1_0, [1, 1, 1] * Dynein.K_OFF_2_0]
        KBT = 4.14195; % pN nm (300 K)
        D0 = 6; % nm
        ALPHA = 0.3;
        BETA = 0.7;
        P_SYN = 0.23;
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
        
        function val = nextAtpOn(obj)
            if obj.PrimaryS < obj.NUM_HYDRO_SITES
                val = obj.PrimaryS + 1;
            elseif obj.S < obj.NUM_SITES
                val = obj.S + 1;
            else
                val = 0;
            end
        end
        function val = kon(obj, atpConc, force)
            nextOn = nextAtpOn(obj);
            kon0 = obj.K_ON_0(nextOn);
            val = kon0 .* atpConc .* obj.loadfactor(force);
        end
        
        function val = nextAtpOff(obj)
            if obj.SecondaryS > 0
                val = obj.NUM_HYDRO_SITES + obj.SecondaryS;
            else
                val = obj.S;
            end
        end
        function val = koff(obj, atpConc)
            nextOff = nextAtpOff(obj);
            koff0 = obj.K_ON_0(nextOff);
            val = koff0 .* atpConc;
        end    
    end % methods
    
    methods (Static)
        function val = stepsize(s)
            val = (Dynein.NUM_SITES - s) * Dynein.STEP_INCREMENT;
        end
        function val = loadfactor(force)
            val = exp((force .* Dynein.D0 ./ (Dynenin.KBT)));
        end
    end % static methods
end % Dynein class