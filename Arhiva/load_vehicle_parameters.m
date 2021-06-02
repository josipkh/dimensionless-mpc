% load vehicle parameters for the chosen model

switch vehicle_params
    case 1
        % parameters IRS 2001
        m = 6.02;  % kg
        Iz = 0.153;  % kgm^2
        Lf = 0.137;  % m
        Lr = 0.222;  % m
        Cf = 40;  % N/rad
        Cr = 52;  % N/rad
        V = 1.98;  % m/s
    case 2
        % parameters full size vehicle
        m = 1670;  % kg
        Iz = 2100;  % kgm^2
        Lf = 0.99;  % m
        Lr = 1.7;  % m
        Cf = 123190;  % N/rad
        Cr = 104190;  % N/rad
        V = 15;  % m/s
end
        
if f ~= 1
    m = m*f;
    Iz = Iz*f^3;
    Lf = Lf*f;
    Lr = Lr*f;
    Cf = Cf*f;
    Cr = Cr*f;
    V = V*sqrt(f);
end
    
