classdef SinghConstants
    properties (Constant)
        NUM_SITES = 4
        NUM_HYDRO_SITES = 1
        STEP_INCREMENT = 8 % nm
        K_ON_1_0 = 4e5 % 1/M 1/s
        K_OFF_1_0 = 10 % 1/s
        K_OFF_2_0 = 250 % 1/s
        K_CAT_0 = 55 % 1/s
        K_ON_0_RATIO = [1, 1, 1/4, 1/6];
        KBT = 4.14195; % pN nm (@ 300 K)
        D0 = 6; % nm
        ALPHA = 0.3;
        P_SYN_0 = 0.23;
        K_TRAP = 7e-3 % pN / nm
    end
    methods (Static)
        function val = restoringForce(x)
            val = SinghConstants.K_TRAP .* x;
        end
    end
end

