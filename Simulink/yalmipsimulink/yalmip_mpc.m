classdef yalmip_mpc < matlab.System & matlab.system.mixin.Propagates
    % An MPC controller implemented using YALMIP
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties(Nontunable)
        nx = 8
        nu = 4
        N = 10
        M = 3
        Ts = 0.05
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)
        Controller
    end

    methods(Access = protected)
        function setupImpl(obj)
            loadStruct = load('SavedController.mat','controller');
            obj.Controller = loadStruct.controller;
        end

        function [predictions,uopt] = stepImpl(obj,Xi,oldUs,xirefs,As,Bs,xihats,deltaf)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            inputs = {Xi,oldUs,xirefs,As,Bs,xihats,deltaf};
            [solutions, errorcode, errortext] = obj.Controller(inputs);
    
            % check the optimization results
            if errorcode ~= 0
                if iscell(errortext)
                    error(errortext{1});
                else
                    errror(errortext);
                end
            end
            predictions = [Xi,solutions{2}(:,2:end)];
            uopt = solutions{1};
        end
        
        function sts = getSampleTimeImpl(obj)
            % Set the sample time
            sts = createSampleTime(obj,'Type','Discrete',...
                                       'SampleTime',obj.Ts,...
                                       'OffsetTime',0);
            %sts = createSampleTime(obj,'Type','Fixed In Minor Step');
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
        function num = getNumInputsImpl(~)
            num = 7;
        end
        function num = getNumOutputsImpl(~)
            num = 2;
        end
        function [dt1,dt2] = getOutputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
        end
        function dt1 = getInputDataTypeImpl(~)
        	dt1 = 'double';
        end
        function [sz1,sz2] = getOutputSizeImpl(obj)
        	sz1 = [obj.nx,obj.N+1];  % predictions
            sz2 = [obj.nu,obj.M];    % optimal inputs
        end
        function sz1 = getInputSizeImpl(~)
        	sz1 = [1,1];
        end
        function cp1 = isInputComplexImpl(~)
        	cp1 = false;
        end
        function [cp1,cp2] = isOutputComplexImpl(~)
        	cp1 = false;
            cp2 = false;
        end
        function fz1 = isInputFixedSizeImpl(~)
        	fz1 = true;
        end
        function [fz1,fz2] = isOutputFixedSizeImpl(~)
        	fz1 = true;
            fz2 = true;
        end
    end
end
