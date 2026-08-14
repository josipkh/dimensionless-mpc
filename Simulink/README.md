# Simulink implementation
This folder contains two Simulink implementations of the dimensionless controller. `SimpleModel` implements the vehicle through a classic ODE, while `VehicleDynamicsBlockset` uses MATLAB's [vehicle dynamics toolbox](https://www.mathworks.com/products/vehicle-dynamics.html).

In both cases, the controller is implemented using YALMIP.