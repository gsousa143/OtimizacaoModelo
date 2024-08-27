classdef OptModel < matlab.System
    % untitled4 Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object.

    % Public, tunable properties
    properties
        limits; np=4; T_otm=15; loose;
    end

    % Pre-computed constants or internal states
    properties (Access = private)
    end

    methods (Access = protected)
        function setupImpl(obj)
            
        end

        function constrains = stepImpl(obj,u,qbot_data,initial_constants,initial_state)
            limits = [
                min(initial_constants+obj.loose,obj.limits(1,:));
                max(initial_constants-obj.loose,obj.limits(2,:))];

            [~, constrains] = de(@(x) costf(x,u,qbot_data,initial_state), ...
                limits, initial_constants, obj.np, Inf, obj.T_otm);
        end

        function resetImpl(obj)
            % Initialize / reset internal properties
        end
        function sts = getSampleTimeImpl(obj)
            sts = createSampleTime(obj,'Type','Discrete','SampleTime',obj.T_otm);
        end
    end

    methods (Access = public)

        function setParams(obj,constantes)
            obj.limits = [
                constantes(1:4)*1.1, 0.07, constantes(6:14)*1.5;
                constantes(1:4)*0.9, 0.00, constantes(6:14)*0.5];

            obj.loose = (obj.limits(1,:)-obj.limits(2,:))*0.01;
        end
    end
end

function cost = costf(x,u,qbot_data,initial_state)
n_data = size(qbot_data,1);
ddmr_data = zeros(n_data,2);
ddmr = DDMR;
ddmr.setup();
ddmr.setState(initial_state);
ddmr.setConstants(x);
for i = 1:n_data
    X = ddmr.step(u(i,:)');
    ddmr_data(i,:) = X(1:2);
end

cost = sqrt(sum((qbot_data-ddmr_data).^2, 'all')/n_data);

end


function [Fbest, gbest, it, x] = de(costf, limites, x0,NP, itMax, tempoMax)

it = 1;
tempoInicial = tic;

D = size(limites,2);

xMax = limites(1,:);
xMin = limites(2,:);


% calcula posicao iniciais e inicializa vetor candidato
x = inicializaPopulacao(xMin,xMax,x0,NP,D);


% calcula o funcional custo
for n = 1:NP
    F(n) = costf(x(n, :));
end
tempoIt = toc(tempoInicial);

Fp = 0.5*ones(NP,1);
cr = 0.9*ones(NP,1);

while(1)
    % determina Fbest e gbest
    [C, I] = min(F);
    Fbest = C;
    gbest = x(I, :);


    %verifica criteiro de parada
    if ( it >= itMax )||( toc(tempoInicial)+tempoIt > tempoMax )
        it
        break;
    end

    v = multacaoRand1(x,NP,Fp);

    v = max(xMin, min(xMax, v));

    u = cruzamento(x,v,NP,D,cr);


    [x,F] = selecao(x,F,u,costf,NP);

    it = it+1;
end



    function x = inicializaPopulacao(xMin,xMax,x0,NP,D)
        for d = 1:D
            x(1:NP,d) = xMin(d) + rand(NP,1)*( xMax(d) - xMin(d) );
        end

        if ~isempty(x0)
            x(1:size(x0,1),:) = max(min(x0,xMax),xMin);
        end
    end


    function v = multacaoRand1(x,NP,Fp)
        v = x(randperm(NP),:) + Fp.*(x(randperm(NP),:) - x(randperm(NP),:) );
    end


    function u = cruzamento(x,v,NP,D,cr)
        i = rand(NP,D)<=cr;

        u = x;
        u(i) = v(i);
    end


    function [x,F] = selecao(x,F,u,costf,NP)
        % calcula o funcional custo dos vetores canditatos
        for n = 1:NP
            Fu(n) = costf(u(n, :));
        end
        % atualiza o vetor de custos
        for n = 1:NP
            if Fu(n) < F(n)
                x(n,:) = u(n,:);
                F(n) = Fu(n);
            end
        end
    end


end


