classdef Dynein
    properties (SetAccess = 'private')
        AtpSites
        Position = 0;
    end % private properties
    
    properties (Constant)
        NUM_SITES = 4
        NUM_HYDRO_SITES = 1
        STEP_INCREMENT = 8
        K_ON_1_0 = 4E5 % 1/M 1/s
        K_OFF_1_0 = 10 % 1/s
        K_OFF_2_0 = 250 % 1/s
        K_CAT_0 = 55 % 1/s
        K_ON_0_RATIO = [1, 1, 1/4, 1/6];
        KBT = 4.14195; % pN nm (@ 300 K)
        D0 = 6; % nm
        ALPHA = 0.3;
        P_SYN_0 = 0.23;
        
        % Dependent Constants
        NUM_SECOND_SITES = Dynein.NUM_SITES - Dynein.NUM_HYDRO_SITES
        
        STEP_SIZE = (Dynein.NUM_SITES + 1 ...
            - colon(1, Dynein.NUM_SITES)) ...
            .* Dynein.STEP_INCREMENT
        
        K_ON_0 = Dynein.K_ON_0_RATIO .* Dynein.K_ON_1_0
        
        K_OFF_0 = [Dynein.K_OFF_1_0, ones(1, Dynein.NUM_SECOND_SITES) ...
            .* Dynein.K_OFF_2_0]
        
        HYDRO_PREFACTOR = [1, ones(1, Dynein.NUM_SECOND_SITES) * 1e-2]
        
        BETA = 1 - Dynein.ALPHA;
    end % constant properties
    
    properties (Dependent)
        S
        PrimaryS
        SecondaryS
    end % dependent properties
    
    methods
        function obj = Dynein
            if Dynein.NUM_SITES <= Dynein.NUM_HYDRO_SITES
                error('NUM_SITES must be greater than NUM_HYDRO_SITES');
            end
            if Dynein.NUM_SITES ~= length(Dynein.K_ON_0_RATIO)
                error('K_ON_0_RATIO must have length equal to NUM_SITES');
            end
            obj.AtpSites = zeros(1, Dynein.NUM_SITES);
        end % Dynein constructor
        
        function val = get.S(obj)
            val = sum(obj.AtpSites);
        end % get.S(obj)
        
        % PrimaryS Get/Set
        function val = get.PrimaryS(obj)
            val = sum(obj.AtpSites(1:Dynein.NUM_HYDRO_SITES));
        end
        function obj = set.PrimaryS(obj, val)
            val = CNSUtils.bound(val, 0, Dynein.NUM_HYDRO_SITES, ...
                'PrimaryS');
            obj.AtpSites(1:end) = false;
            obj.AtpSites(1:val) = true;
        end
        
        % SecondaryS Get/Set
        function val = get.SecondaryS(obj)
            val = sum(obj.AtpSites((Dynein.NUM_HYDRO_SITES + 1):end));
        end
        function obj = set.SecondaryS(obj, val)
            val = CNSUtils.bound(val, 0, Dynein.NUM_SECOND_SITES, ...
                'SecondaryS');
            iBegin = Dynein.NUM_HYDRO_SITES + 1;
            obj.AtpSites(iBegin:end) = false;
            obj.AtpSites(iBegin:val) = true;
        end
        
        % AtpSites Setter
        function obj = set.AtpSites(obj, val)
            if length(val) ~= obj.NUM_SITES
                error('AtpSites must be a vector of length %d.', ...
                    Dynein.NUM_SITES);
            end
            logicalVal = logical(val);
            if ~all(val == logicalVal)
                warning('AtpSites muct be a logical vector.');
                val = logicalVal;
            end
            obj.AtpSites = val;
        end % set.AtpSites(obj, val)
        
        % Add and remove order
        function val = nextAtpOn(obj)
            val = find(~obj.AtpSites, 1);
        end
        function val = nextAtpOff(obj)
            val = find(obj.AtpSites, 1, 'last');
        end
        function val = nextAtpHydro(obj)
            primarySites = obj.AtpSites(1:Dynein.NUM_HYDRO_SITES);
            val = find(primarySites, 1, 'last');
        end
        
        % Adders and removers
        function obj = bindAtp(obj)
            nextOn = obj.nextAtpOn;
            if isempty(nextOn)
                warning('Sites full, cannot bind another ATP.');
            end
            obj.AtpSites(nextOn) = true;
        end
        function obj = unbindAtp(obj)
            nextOff = obj.nextAtpOff;
            if isempty(nextOff)
                warning('Sites empty, cannot unbind another ATP.');
            end
            obj.AtpSites(nextOff) = false;
        end
        function obj = hydrolyzeAtp(obj)
            nextHydro = obj.nextAtpHydro;
            if isempty(nextHydro)
                warning(['Primary sites empty, cannot '
                    'hydrolyze another ATP.']);
            end
            obj.AtpSites(nextHydro) = false;
        end
        
        % Rates Constants
        function val = kon(obj, atpConc, force)
            nextOn = nextAtpOn(obj);
            if isempty(nextOn)
                val = 0;
                return
            end
            kon0 = Dynein.K_ON_0(nextOn);
            s = obj.S;
            if s <= Dynein.NUM_HYDRO_SITES
                loadfactor = 1;
            else
                loadfactor = Dynein.loadfactor(force);
            end
            val = kon0 .* atpConc .* loadfactor;
        end
        function val = koff(obj, atpConc)
            nextOff = nextAtpOff(obj);
            if isempty(nextOff)
                val = 0;
                return
            end
            koff0 = Dynein.K_ON_0(nextOff);
            val = koff0 .* atpConc;
        end
        function val = kcat(obj, force)
            if isempty(obj.nextAtpHydro)
                val = 0;
                return;
            end
            s = obj.S;
            prefactor = Dynein.HYDRO_PREFACTOR(s);
            expTerm = exp(-Dynein.ALPHA .* force .* Dynein.stepsize(s) ...
                ./ Dynein.KBT);
            val = prefactor .* Dynein.K_CAT_0 .* expTerm;
        end
        function val = getRateArray(obj, atpConc, force)
            val = [obj.kon(atpConc, force), ...
                   obj.koff(atpConc), ...
                   obj.kcat(force)];
        end
        
        % stepping
        function val = probabilityReverse(obj, force)
            expTerm = exp(Dynein.BETA .* force .* obj.currentStepSize ...
                          ./ Dynein.KBT);
            val = Dynein.P_SYN_0 .* expTerm;
            val(val > 1) = 1;
        end
        function attemptStep(obj)
            
        end
        % state dependent calculations
        function val = currentStepSize(obj)
            s = obj.S;
            val = Dynein.stepsize(s);
        end

    end % methods
    
    methods (Static)
        function val = stepsize(s)
            if s == 0, val = 0; return; end
            val = Dynein.STEP_SIZE(s);
        end
        function val = loadfactor(force)
            val = exp((force .* Dynein.D0 ./ (Dynein.KBT)));
        end
        
    end % static methods
end % Dynein class