% load vehicle parameters for the chosen model

switch vehicle_params
    case 1
        % parameters IRS 2001
        m = 6.02;  % kg
        Iz = 0.153;  % kgm^2
        Lf = 0.137;  % m
        Lr = 0.222;  % m
        Cf = 43*2;  % N/rad (2 tires)
        Cr = 43*2;  % N/rad (2 tires)
        V = 1.98;  % m/s
    case 2
        % parameters full size (Carmaker demo compact)
        m = 1194;  % kg
        Iz = 1645.687;  % kgm^2
        Lf = 1.01;  % m
        Lr = 1.58;  % m
        Cf = 67043*2;  % N/rad (2 tires)
        Cr = 67043*2;  % N/rad (2 tires)
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
