classdef Dynein
% DYNEIN A class representing a dynein motor.    
    
    
    
    properties (SetAccess = 'private')
        AtpSites
        Position = 0;
    end % private properties
    
    
    
    properties (Constant)        
        % From Paper
        NUM_SITES =  SinghConstants.NUM_SITES
        NUM_HYDRO_SITES = SinghConstants.NUM_HYDRO_SITES
        STEP_INCREMENT = SinghConstants.STEP_INCREMENT
        K_ON_1_0 = SinghConstants.K_ON_1_0
        K_OFF_1_0 = SinghConstants.K_OFF_1_0
        K_OFF_2_0 = SinghConstants.K_OFF_2_0
        K_CAT_0 = SinghConstants.K_CAT_0
        K_ON_0_RATIO = SinghConstants.K_ON_0_RATIO
        KBT = SinghConstants.KBT
        D0 = SinghConstants.D0
        ALPHA = SinghConstants.ALPHA
        P_SYN_0 = SinghConstants.P_SYN_0
        
        % Dependent Constants
        NUM_SECOND_SITES = Dynein.NUM_SITES - Dynein.NUM_HYDRO_SITES
        
        STEP_SIZE = (Dynein.NUM_SITES + 1 ...
            - colon(1, Dynein.NUM_SITES)) ...
            .* Dynein.STEP_INCREMENT
        
        K_ON_0 = Dynein.K_ON_0_RATIO .* Dynein.K_ON_1_0
        
        K_OFF_0 = [Dynein.K_OFF_1_0, ones(1, Dynein.NUM_SECOND_SITES) ...
            .* Dynein.K_OFF_2_0]
        
        HYDRO_PREFACTOR = [1e-2, ones(1, Dynein.NUM_SECOND_SITES)]
        
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
                return;
            end
            kon0 = Dynein.K_ON_0(nextOn);
            if nextOn > Dynein.NUM_HYDRO_SITES
                loadfactor = Dynein.loadfactor(force);                
            else
                loadfactor = 1;
            end
            val = kon0 .* atpConc .* loadfactor;
        end
        function val = koff(obj)
            nextOff = nextAtpOff(obj);
            if isempty(nextOff)
                val = 0;
                return;
            end
            koff0 = Dynein.K_OFF_0(nextOff);
            val = koff0;
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
        function val = getRates(obj, atpConc, force)
            val = [obj.kon(atpConc, force), ...
                obj.koff, ...
                obj.kcat(force)];
        end
        
        % stepping
        function val = probabilityReverse(obj, force)
            s = obj.S;
            expTerm = exp(Dynein.BETA .* force .* Dynein.stepsize(s) ...
                ./ Dynein.KBT);
            val = Dynein.P_SYN_0 .* expTerm;
            val = CNSUtils.bound(val, 0, 1);
        end
        function obj = step(obj)
            obj.Position = obj.Position + Dynein.stepsize(obj.S);
            obj = obj.hydrolyzeAtp;
        end
        function obj = attemptStep(obj, force)
            if ~CNSUtils.randchoose(obj.probabilityReverse(force))
                obj = obj.step;
            end
        end
        
        function [obj, dt] = update(obj, atpConc, force)
            sz = size(obj);
            nObj = prod(sz);
            dt = zeros(sz);
            % FIXME - Functionalize this
            if isscalar(atpConc)
                atpConc = atpConc .* ones(sz);
            end
            if isscalar(force)
                force = force .* ones(sz);
            end
            for i = 1:nObj
                currentObj = obj(i);
                rates = currentObj.getRates(atpConc(i), force(i));
                dt(i) = CNSUtils.mcTimeStep(rates); 
                choice = CNSUtils.randchoose(rates); 
                switch choice
                    case 1
                        obj(i) = currentObj.bindAtp;
                    case 2
                        obj(i) = currentObj.unbindAtp;
                    case 3
                        obj(i) = currentObj.attemptStep(force(i));            
                end % switch
            end % loop through objects
        end % update function
        
        function disp(obj)
            siteFormat = [repmat('%d - ', 1, Dynein.NUM_SITES - 1), '%d'];
            siteStr = sprintf(siteFormat, obj.AtpSites);
            fprintf(['\tDynein with ATP Sites: %s | ', ...
                     'at Position: %5.1f\n'], ... 
                    siteStr, ...
                    obj.Position);
        end
    end % methods
    
    methods (Static)
        function val = stepsize(s)
            if s == 0, val = 0; return; end
            val = Dynein.STEP_SIZE(s);
        end
        
        function val = loadfactor(force)
            val = exp((force .* Dynein.D0 ./ Dynein.KBT));
        end
    end % static methods
end % Dynein class