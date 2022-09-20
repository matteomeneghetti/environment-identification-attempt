classdef FilterNthOrder < handle
    %FILTERNTHORDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Access = private )
        h
        alfa 
        beta 
        beta_hat 
        X 
        X_old
        a
        b
        order
    end
    
    methods
        function obj = FilterNthOrder(order, a, b, h)
            obj.order = order;
            obj.a = a;
            obj.b = b;
            obj.h = h;
        end
        
        function y = process(obj, u, t)
            %% Initialise.
            if isempty(obj.alfa)
                obj.alfa = zeros(1, obj.order + 1);
                obj.beta = zeros(1, obj.order + 1);

                % Scattolini notation:
                for i = 1 : (obj.order + 1)
                    obj.alfa(i) = obj.a(i) / obj.a(1);
                    obj.beta(i) = obj.b(i) / obj.a(1);
                end

                % DENOMINATOR: alfa(n) * s ^ n + ... + alfa(1)
                obj.alfa = fliplr(obj.alfa);
                % NUMERATOR: beta(n) * s ^ n + ... + beta(1)
                obj.beta = fliplr(obj.beta);

                % Symplectic parameters
                obj.beta_hat = zeros(1, obj.order + 1);

                % b_hat_{n - 1}
                % Es.: order = 2 -> b_hat_2 = B_HAT(3) = beta(3) = beta_2
                obj.beta_hat(obj.order + 1) = obj.beta(obj.order + 1);
                % All the others.
                for i = 1 : obj.order
                    obj.beta_hat(i) = obj.beta(i) - obj.alfa(i) * obj.beta_hat(obj.order + 1);
                end
            end


            %% Filter
            if isempty(obj.X)
                obj.X = zeros(1, obj.order);
                obj.X_old = zeros(1, obj.order);
            end

            F = zeros(1, obj.order);

            for i = 1 : obj.order
                F(obj.order) = F(obj.order) - obj.alfa(i) * obj.X(i);
            end

            F(obj.order) = F(obj.order) + u;
            obj.X(obj.order) = obj.X_old(obj.order) + obj.h * F(obj.order);
            obj.X_old(obj.order) = obj.X(obj.order);

            y = obj.beta_hat(obj.order) * obj.X(obj.order) + obj.beta_hat(obj.order + 1) * u;

            for i = (obj.order - 1) : -1 : 1
                F(i) = obj.X(i + 1);
                obj.X(i) = obj.X_old(i) + obj.h * F(i);
                obj.X_old(i) = obj.X(i);

                y = y + obj.beta_hat(i) * obj.X(i);
            end
        end
    end
end

