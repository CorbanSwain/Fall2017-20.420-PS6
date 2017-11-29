classdef Vesicle < handle
    
    
    properties
        
        Position
        Dynein
        
    end
    
    properties (Dependent)
        isbound
    end
    
    methods
        
        function objs = Vesicle(startPos)
            if nargin > 0
                nVes = size(startPos, 1);
                objs(1, nVes) = Vesicle;
                for i = 1:nVes
                    objs(i).Position = startPos(i, :);                    
                end
            end
        end
        
        function val = get.isbound(obj)
            val = ~isempty(obj.Dynein);
        end
        
        function objs = update(objs, dt, time, dyneins)
            
            for iObj = 1:numel(objs)
                obj = objs(iObj);
                if obj.isbound
                    obj.Position = obj.Dynein.interpPos(time) ...
                        + [-ModConstants.MAX_DIST 0];
                    poff = ModConstants.K_OFF_MEM_DYN .* dt;
                    % May need with different dt
                    % poff = CNSUtils.bound(poff,0,1);
                    if CNSUtils.randchoose(poff)
%                         fprintf('\nUN-BINDING,        Ves-%2d\n', iObj);
                        obj.Dynein.Vesicle = [];
                        obj.Dynein = [];
                    end
                else % not bound
                    obj.Position = rand_walk(obj.Position, ...
                        ModConstants.D, ...
                        dt, ModConstants.SIM_X_LIM,...
                        ModConstants.SIM_Y_LIM);
                    if ~isempty(dyneins)
                        distances = obj.distTo(time, dyneins);
                        pons = ModConstants.K_ON_EFF(distances) .* dt;
                        if any(pons)
                            pSomething = CNSUtils.bound(sum(pons), 0, 1);
                            if CNSUtils.randchoose(pSomething)
                                choice = CNSUtils.randchoose(pons);
%                                 fprintf('\nBINDING, Dyn-%2d to Ves-%2d\n', ...
%                                     choice, iObj);
                                obj.Dynein = dyneins(choice);
                                obj.Dynein.Vesicle = obj;
                            end
                        end
                    end
                end
                objs(iObj) = obj;
            end % main for loop
        end % update funciton
        
        function val = distTo(obj, time, dyneins)
            p = obj.Position;
            sz = size(dyneins);
            nDyn = prod(sz);
            val = zeros(sz);
            for i = 1:nDyn
                dyn = dyneins(i);
                if dyn.HasVesicle
                    % ^ this should technically be handled elsewhere ...
                    val(i) = ModConstants.MAX_DIST;
                else
                    dynPos = dyn.interpPos(time);
                    delta = p -  dynPos;
                    if any(delta >= ModConstants.MAX_DIST)
                        val(i) = ModConstants.MAX_DIST;
                    else
                        val(i) = sqrt(sumsqr(delta));
                    end
                end
            end
        end % distTo function
        
    end
    
end