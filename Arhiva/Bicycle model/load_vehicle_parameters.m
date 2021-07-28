% load vehicle parameters for the chosen model

switch vehicle_params
    case 1
        % parameters lab size
        load ../../ParamsLab
        V = 60/3.6/10 - 0.1;  % m/s
    case 2
        % parameters full size
        load ../../ParamsFull
        V = 60/3.6;  % m/s
end

m = VEHICLE.MASS;
Iz = VEHICLE.INERTIA_Z;
Lf = VEHICLE.LF;
Lr = VEHICLE.LR;
Cf = VEHICLE.CORNERING_STIFF;  % single tire
Cr = VEHICLE.CORNERING_STIFF;
        
if f ~= 1
    m = m*f;
    Iz = Iz*f^3;
    Lf = Lf*f;
    Lr = Lr*f;
    Cf = Cf*f;
    Cr = Cr*f;
    V = V*sqrt(f);
end
