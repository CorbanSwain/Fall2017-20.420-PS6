classdef ModConstants
    properties (Constant)
        ETA = 6E3 % pN s / nm^2
        R_VES = 65 % nm
        L_DYN = 35 % nm
        MEM_AREA_DENSITY = 4.23E-24 % mol membrane / nm^2
        K_D_MEM_DYN = 6E3 % pM
        
        
        ZETA = 6 * pi * ModConstants.ETA * ModConstants.R_VES % pN s / nm
        D = SinghConstants.KBT / ModConstants.ZETA % nm^2 / s
        VD = 4 / 3 * pi * ModConstants.L_DYN .^ 3
        MAX_DIST = ModConstants.R_VES + ModConstants.L_DYN
        MIN_DIST = abs(ModConstants.R_VES - ModConstants.L_DYN)
    end
    
    methods (Static)
        function f = drag_force(v)
            f = ModConstants.ZETA .* v; % v in nm / s, F in pN
        end
        
        function val = SA_V(d, check)
            % SA_V Surface area of the vessicle within dynein's effective
            % volume.
            %
            % units - nanometers squared
            %
            % d: the distance from dynein's location to the center of the
            % vessicle
            
            
            switch nargin
                case 1
                    check = true;
                case 2
                otherwise
                    error('Unexpected number of arguments');
            end
            if check
                d = abs(d);
                sz = size(d);
                val = zeros(sz);
                nD = prod(sz);
                for i = 1:nD
                    if d(i) <= ModConstants.MIN_DIST
                        val(i) = 0;
                    elseif d(i) >= ModConstants.MAX_DIST
                        val(i) = 0;
                    else
                        % Find this derivation by looking through the
                        % equation helper
                        val(i) = -(65.*pi.*(d(i).^2 - 130.*d(i) + 3000))./d(i);
                    end
                end
            else
                val = -(65.*pi.*(d.^2 - 130.*d + 3000))./d;
            end
        end
        
        function val = VD_EFF(d, check)
            % VD_EFF effective volume that dynein can move in.
            %
            % units - nanometers cubed
            %
            % this value collapses to zero as a vessicle moves closer to
            % mimic the effect of hightened binding probability as a
            % vessicle comes into close contact with the dynein
            
            
            switch nargin
                case 1
                    check = true;
                case 2
                otherwise
                    error('Unexpected number of arguments');
            end
            if check
                d = abs(d);
                sz = size(d);
                val = zeros(sz);
                nD = prod(sz);
                for i = 1:nD
                    if d(i) <= ModConstants.MIN_DIST
                        val(i) = 0;
                    elseif d(i) >= ModConstants.MAX_DIST
                        val(i) = ModConstants.VD;
                    else
                        % Find this derivation by looking through the
                        % equation helper
                        val(i) = -(pi.*(d(i) - 30).^2 ...
                            .*(d(i).^2 + 60.*d(i) - 30000)) ...
                            ./(12.*d(i));
                    end
                end
            else
                val = -(pi.*(d - 30).^2 ...
                    .*(d.^2 + 60.*d - 30000)) ...
                    ./(12.*d);
            end
        end
        
        function val = memConc(d)
            % Effective concentration of membrane in proximity to dynein
            %
            % units - molar
            
            sz = size(d);
            val = zeros(sz);
            nD = prod(sz);
            for i = 1:nD
                if d(i) >= ModConstants.MAX_DIST
                    val(i) = 0;
                elseif d(i) <= ModConstants.MIN_DIST
                    val(i) = Inf;
                else
                    val(i) = ModConstants.MEM_AREA_DENSITY ...
                        .* ModConstants.SA_V(d(i), false) ... 
                        ./ ModConstants.VD_EFF(d(i), false) * 1e24;
                end
            end
        end
        
        function equationHelper
            syms d rd rv;
            hv = (rd - rv + d) .* (rd + rv - d) ./ (2 .* d);
            a = 1 ./ (2 .* d) .* sqrt(4 .* d.^2 .* rd.^2 ...
                - (d.^2 - rv.^2 + rd.^2).^2);
            SA_v = pi .* (hv .^ 2 + a .^ 2);
            SA_v = simplify(subs(SA_v,[rd rv], ...
                [ModConstants.L_DYN, ...
                ModConstants.R_VES]));
            display(SA_v);
            
            Vint = pi ./ (12 * d) .* (rd + rv - d) .^ 2 ...
                .* (d.^2 + 2*d*rv - 3*rv.^2 ...
                + 2.*d.*rd + 6.*rv.*rd - 3.*rd.^2);
            Vd_eff = ModConstants.VD - Vint;
            Vd_eff = simplify(subs(Vd_eff,[rd rv], ...
                [ModConstants.L_DYN, ...
                ModConstants.R_VES]));
            display(Vd_eff);
            
        end
    end
end