classdef Vesicle < handle
    
    
    properties
        
        Position
        Time = 0
        Dynein
        
    end
    
    
    methods
        
        function obj = Vessicle(startPos)
            obj.Position = startPos;
        end
        
        function val = isbound(obj)
            val = ~isempty(obj.Dynein);
        end
        
        function objs = update(objs, dt, time, dyneins)
            
            for iObj = 1:numel(objs)
                obj = objs(iObj);
                obj.Time = obj.Time + dt;
                if obj.isbound
                    obj.Position = obj.Dynein.interpPos(time) ...
                        + [-ModConstants.MAX_DIST 0];
                    poff = ModConstants.K_OFF_MEM_DYN .* dt;
                    % May need with different dt
                    % poff = CNSUtils.bound(poff,0,1);
                    if CNSUtils.randchoose(poff)
                        obj.Dynein = [];
                    end
                else % not bound
                    obj.Position = rand_walk(obj.Position, ...
                        ModConstants.D, ...
                        dt, ModConstants.SIM_X_LIM,...
                        ModConstants.SIM_Y_LIM);
                    distances = obj.distTo(time, dyneins);
                end
                objs(iObj) = obj;
            end % main for loop
        end
        
        function val = distTo(obj, time, dyneins)
            p = obj.Position;
            for i = 1:numel(dyneins)
                dyn = dyneins(i);
                dynPos = dyn.interpPos(time);
                deltap = dyn.
            end
        end
        
    end
    
end