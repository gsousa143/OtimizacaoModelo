classdef DDMR < matlab.System %#codegen
    % untitled Add summary here
    %
    % MATLAB Object System Modelado para atualizar o estado de um DDMR
    % utilizando o modelo unificado de acordo com as constatens setadas.

    % Public, tunable properties
    properties
        R = 0; L = 0; alpha_s= 0 ; alpha_k=0; f_s=0; f_k=0; V=zeros(2,2); A_tau=zeros(2,2); B_tau=zeros(2,2); M_i=zeros(2,2); k_i=0; k_p=0; T = 0.01;
    end

    % Pre-computed constants or internal states
    properties (Access = private)
        state;
    end

    methods (Access = protected)
        function setupImpl(obj)
            obj.state = zeros(7,1);
            % % Perform one-time calculations, such as computing constants
        end
        

        function state = stepImpl(obj,u)

            % Recebe um sinal de controle e atualiza o estado do ddmr
            [~, states] = ode23(@(t,X) ...
                unifiedModel(t, X, u, obj.R, obj.L, obj.alpha_s, obj.alpha_k, obj.f_s, obj.f_k, obj.V, obj.A_tau, obj.B_tau, obj.M_i, obj.k_i, obj.k_p), ...
                [0,obj.T], obj.state);

            state = states(end,:);
            obj.state = state';
        end

        function resetImpl(obj)
            % Initialize / reset internal properties
        end


        function [o1] = getOutputSizeImpl(obj)
            % Return size for output port to be same as input port
            o1 = [1 7];

        end

        function [o1] = getOutputDataTypeImpl(obj)
            % Return data type for output port to be same as input port
            o1 = "double";
 
        end

        function [o1] = isOutputComplexImpl(~)
            % Return output port complexity to be real
            o1 = false;

        end

        function [o1] = isOutputFixedSizeImpl(~)
            % Return true for each output port with fixed size
            o1 = true;

        end

        function sts = getSampleTimeImpl(obj)
            sts = createSampleTime(obj,'Type','Discrete','SampleTime',obj.T);
        end
    end

    methods(Access = public)
        function state = getState(obj)
            state = obj.state;
        end

        function setState(obj,state)
            obj.state = state;
        end


        function setConstants(obj, constants)
            obj.R = constants(1); % Raio das rodas do robô (m)
            obj.L = constants(2); % Comprimento do semieixo das rodas do robô (m)
            m_c = constants(3); % Massa da plataforma do robô
            m_w = constants(4); % Massa da roda do robô
            d = constants(5);

            N = constants(6); %Relacao de engrenagens do motor CC
            b = constants(7); % Coeficiente de atrito viscoso do motor CC
            K_ce = constants(8); % Constante de força contra-eletromotriz do motor CC
            K_t = constants(9); % Constante de torque do motor CC
            R_a = constants(10); % Resistencia de armadura totau do robô

            f_s = constants(11); % Coeficiente de atrito estático
            f_k = constants(12); % Coeficiente de atrito cinético

            obj.alpha_s = constants(13); % Constante de saturação do atrito estático
            obj.alpha_k = constants(14); % constants de saturação do atrito cinético

            obj.k_i = 1; % Ganho Ki
            obj.k_p = 2.5; % Ganho Kp


            g = 9.81; %Aceleracao da gravidade

            m = m_c + 2*m_w;
            I_c = m_c*obj.L^2/2;
            I_w = m_w*obj.R^2/2;
            I_t = I_c + 2*m_w*obj.L^2 + 2*I_w + m_c*d^2;

            obj.f_s = f_s*m*g;
            obj.f_k = f_k*m*g;

            obj.V = obj.R^2*d*m_c/(4*obj.L^2)*[
                0, 1;
                -1,0];

            obj.A_tau = (b - (K_ce*K_t)/R_a)/N*eye(2,2);

            obj.B_tau = K_t/R_a*eye(2,2);

            obj.M_i = -inv([
                I_w + obj.R^2*(m*obj.L^2 + I_t)/(4*obj.L^2), obj.R^2*(m*obj.L^2 - I_t)/(4*obj.L^2);
                obj.R^2*(m*obj.L^2 - I_t)/(4*obj.L^2), I_w + obj.R^2*(m*obj.L^2 + I_t)/(4*obj.L^2)]);
        end

    end
end

%Helper Functions
function dstate_dt = unifiedModel(~, state, u, R, L, alpha_s, alpha_k, f_s, f_k, V, A_tau, B_tau, M_i, k_i, k_p)
A_cin = R/2*[
    cos(state(3)), cos(state(3));
    sin(state(3)), sin(state(3))];


tau_F = [
    f_s*tanh(alpha_s*state(4)) - f_k*tanh(alpha_k*state(4));
    f_s*tanh(alpha_s*state(5)) - f_k*tanh(alpha_k*state(5))];

eta = [
    state(4);
    state(5)];


epi = [
    state(6);
    state(7)];

Y_pi = min(12, max( k_i*epi + k_p*(u-eta), -12) );


dtheta_dt = R*(state(4)-state(5))/(2*L);

dstate_dt = [
    A_cin*eta;
    dtheta_dt;
    M_i*((A_tau+(dtheta_dt*V))*eta - B_tau*Y_pi + tau_F);
    u-eta];
end
