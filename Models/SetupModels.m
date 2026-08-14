function [] = SetupModels(VEHICLE, CONST)
% generates MATLAB functions for faster evaluation of vehicle dynamics and
% its Jacobians
if ~exist("./Codegen", "dir")
    mkdir("Codegen")  % ignored in git, used for storing generated functions
end

if ~exist("./Codegen/PredictionModel.m", "file") || ~exist("./Codegen/PredictionModelJacobian.m", "file")
    disp("Dimensional prediction model not found, recreating...")
    FormulatePredictionModel(VEHICLE, CONST);
end

% if ~exist("./Codegen/PredictionModelPi.m", "file") || ~exist("./Codegen/PredictionModelPiJacobian.m", "file")
%     disp("Dimensionless prediction model not found, recreating...")
%     FormulatePredictionModelPi(VEHICLE, CONST);
% end

if ~exist("./Codegen/SimulationModel.m", "file")
    disp("Simulation model not found, recreating...")
    FormulateSimulationModel(VEHICLE, CONST);
end

end