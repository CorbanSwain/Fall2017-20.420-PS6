classdef ModDynein < handle
% DYNEIN A class representing a dynein motor.    
    
    properties
        Vesicle
    end
    
    properties (SetAccess = 'private')
        AtpSites
        Position
        Time = 0
        PrevTime
        PrevPos

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
        NUM_SECOND_SITES = ModDynein.NUM_SITES - ModDynein.NUM_HYDRO_SITES
        
        STEP_SIZE = (ModDynein.NUM_SITES + 1 ...
            - colon(1, ModDynein.NUM_SITES)) ...
            .* ModDynein.STEP_INCREMENT
        
        K_ON_0 = ModDynein.K_ON_0_RATIO .* ModDynein.K_ON_1_0
        
        K_OFF_0 = [ModDynein.K_OFF_1_0, ones(1, ModDynein.NUM_SECOND_SITES) ...
            .* ModDynein.K_OFF_2_0]
        
        HYDRO_PREFACTOR = [1e-2, ones(1, ModDynein.NUM_SECOND_SITES)]
        
        BETA = 1 - ModDynein.ALPHA;
        
    end % constant properties
    
    
    
    properties (Dependent)
        S
        PrimaryS
        SecondaryS
        HasVesicle
        Velocity
    end % dependent properties
    
    
    
    methods
        function objs = ModDynein(startPos)
            % CONSTANT VALIDITY CHECKING
            if ModDynein.NUM_SITES <= ModDynein.NUM_HYDRO_SITES
                error('NUM_SITES must be greater than NUM_HYDRO_SITES');
            end
            if ModDynein.NUM_SITES ~= length(ModDynein.K_ON_0_RATIO)
                error('K_ON_0_RATIO must have length equal to NUM_SITES');
            end
 
            if nargin > 0
                % INITIALIZATION           
                nDyns = size(startPos, 1);
                objs(1, nDyns) = ModDynein;
                for i = 1:nDyns
                    objs(i).Position = startPos(i, :);
                    objs(i).AtpSites = zeros(1, ModDynein.NUM_SITES);
                end
            end
        end % ModDynein constructor
        
        function val = get.S(obj)
            val = sum(obj.AtpSites);
        end % get.S(obj)
        
        % PrimaryS Get/Set
        function val = get.PrimaryS(obj)
            val = sum(obj.AtpSites(1:ModDynein.NUM_HYDRO_SITES));
        end
        function obj = set.PrimaryS(obj, val)
            val = CNSUtils.bound(val, 0, ModDynein.NUM_HYDRO_SITES, ...
                'PrimaryS');
            obj.AtpSites(1:end) = false;
            obj.AtpSites(1:val) = true;
        end
        
        % SecondaryS Get/Set
        function val = get.SecondaryS(obj)
            val = sum(obj.AtpSites((ModDynein.NUM_HYDRO_SITES + 1):end));
        end
        function obj = set.SecondaryS(obj, val)
            val = CNSUtils.bound(val, 0, ModDynein.NUM_SECOND_SITES, ...
                'SecondaryS');
            iBegin = ModDynein.NUM_HYDRO_SITES + 1;
            obj.AtpSites(iBegin:end) = false;
            obj.AtpSites(iBegin:val) = true;
        end
        
        % AtpSites Setter
        function obj = set.AtpSites(obj, val)
            if length(val) ~= obj.NUM_SITES
                error('AtpSites must be a vector of length %d.', ...
                    ModDynein.NUM_SITES);
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
            primarySites = obj.AtpSites(1:ModDynein.NUM_HYDRO_SITES);
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
            kon0 = ModDynein.K_ON_0(nextOn);
            if nextOn > ModDynein.NUM_HYDRO_SITES
                loadfactor = ModDynein.loadfactor(force);                
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
            koff0 = ModDynein.K_OFF_0(nextOff);
            val = koff0;
        end
        function val = kcat(obj, force)
            if isempty(obj.nextAtpHydro)
                val = 0;
                return;
            end
            s = obj.S;
            prefactor = ModDynein.HYDRO_PREFACTOR(s);
            expTerm = exp(-ModDynein.ALPHA .* force .* ModDynein.stepsize(s) ...
                ./ ModDynein.KBT);
            val = prefactor .* ModDynein.K_CAT_0 .* expTerm;
        end
        function val = getRates(obj, atpConc, force)
            val = [obj.kon(atpConc, force), ...
                obj.koff, ...
                obj.kcat(force)];
        end
        
        % stepping
        function val = probabilityReverse(obj, force)
            s = obj.S;
            expTerm = exp(ModDynein.BETA .* force .* ModDynein.stepsize(s) ...
                ./ ModDynein.KBT);
            val = ModDynein.P_SYN_0 .* expTerm;
            val = CNSUtils.bound(val, 0, 1);
        end
        function obj = step(obj)
            obj.Position(1) = obj.Position(1) + ModDynein.stepsize(obj.S);
            obj = obj.hydrolyzeAtp;
        end
        function obj = attemptStep(obj, force)
            if ~CNSUtils.randchoose(obj.probabilityReverse(force))
                obj = obj.step;
            end
        end
        
        function obj = update(obj, atpConc)
            nObj = numel(obj);
            for i = 1:nObj
                if obj(i).HasVesicle
                    force = ModConstants.drag_force(obj(i).Velocity);
                else
                    force = 0;
                end
                
                rates = obj(i).getRates(atpConc, force);
                dt = CNSUtils.mcTimeStep(rates);
                obj(i).PrevTime = obj(i).Time;
                obj(i).Time = obj(i).PrevTime + dt;                                
                choice = CNSUtils.randchoose(rates); 
                
                obj(i).PrevPos = obj.Position;
                switch choice
                    case 1
                        obj(i) = obj(i).bindAtp;                        
                    case 2
                        obj(i) = obj(i).unbindAtp;                        
                    case 3
                        obj(i) = obj(i).attemptStep(force);  
%                         fprintf('STEPPING Dynein %d\n', i);
                end % switch
                
            end % loop through objects
        end % update function
        
        function val = interpPos(obj, tq)
            if isempty(obj.PrevPos)
               val = obj.Position;
               return
            end
            x = [obj.PrevPos(1), obj.Position(1)];
            if x(1) == x(2)
                val = obj.Position;
                return
            end
            t = [obj.PrevTime, obj.Time];
            val = [interp1(t, x, tq), obj.Position(2)];
            if isnan(val(1))
                error(['Querried time is outside of ', ...
                       'the interpolation range.']);
            end
        end
        
        function val = get.HasVesicle(obj)
            val = ~isempty(obj.Vesicle);
        end
        
        function val = get.Velocity(obj)
            if isempty(obj.PrevPos)
               val = 0;
               return
            end
            dx = obj.Position(1) - obj.PrevPos(1);
            if dx == 0
                val = 0;
            else
                val = dx ./ (obj.Time - obj.PrevTime);
            end
        end
        
%         function disp(obj)
%             siteFormat = [repmat('%d - ', 1, ModDynein.NUM_SITES - 1), '%d'];
%             siteStr = sprintf(siteFormat, obj.AtpSites);
%             fprintf(['\tModDynein with ATP Sites: %s | ', ...
%                      'at Position: %5.1f\n'], ... 
%                     siteStr, ...
%                     obj.Position);
%         end

    end % methods
    
    methods (Static)
        function val = stepsize(s)
            if s == 0, val = 0; return; end
            val = ModDynein.STEP_SIZE(s);
        end
        
        function val = loadfactor(force)
            val = exp((force .* ModDynein.D0 ./ ModDynein.KBT));
        end
    end % static methods
end % ModDynein class