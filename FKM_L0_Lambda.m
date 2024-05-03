
function [U, P, J, XieBeni] = FKM_L0_Lambda(X, C, m, lambda, conv, Max_iter, stand, alpha)

%
% Input
% X:        data matrix
% C:        number of clusters
% m:        parameter of fuzziness (e.g. = 2)
% lambda:   parameter of the L0 regularization term (>= 0)
% conv:     convergence criterion (e.g. = 1e-6)
% Max_iter: maximum number of iterations (e.g. = 1e+3)
% stand:    standardization (if = 1, data are standardized)
%
% Output
% U:        membership degree matrix
% P:        prototype matrix
% J:        loss function vector (last element = value at convergence) 
%

[N, T] = size(X); 

    %% standardization
    if stand == 1
	    % Jm = eye(n) - (1/n)*ones(n);
        Jm = eye(N) - (1/N)*ones(N);
	    % X = Jm*X/diag(std(X,1));
        sd = diag(std(X,1));
	    X = Jm*X/sd;
    end

    %% Initialization
	
	P = zeros(C,T);
	D = zeros(N,C);
    Y = rand(N,C);
    U = zeros(N,C);
    for j = 1:N
        U(j,:) = Y(j,:)/sum(Y(j,:));
    end
    
	%% Optimization

	iter = 0;
	J = zeros(Max_iter, 1);
    while iter < Max_iter
		iter = iter + 1;
        % Update of the prototypes 
        for i = 1:C
			P(i,:) = ((U(:,i).^m)'*X)/sum((U(:,i).^m));
        end
        % Compute units-centroids distancies
        epsilon = 1e-10;
		for j = 1:N
			for i = 1:C
                D(j,i) = max(sum((X(j,:)-P(i,:)).^2), epsilon);
            end
        end

        % Update of the membership degrees
		SUM = sum((1./D).^(1/(m-1)),2);
        for j = 1:N
            for i = 1:C
                U(j,i) = (1/(D(j,i))).^(1/(m-1))/SUM(j);
            end

            % Check over U matrix
            checkrj = 1;
            Uopt = U(j,:);
            U0j = U(j,:); 
            p0j = (U0j > 0);
            while checkrj == 1
                Uaux = U0j; 
                Uaux(not(p0j)) = 2;
                [~,mj] = min(Uaux);
                p0j(mj) = 0;
                U0j(mj) = 0;
                SUM0 = sum(p0j.*(1./D(j,:)).^(1/(m-1)),2);
                for i = 1:C 
                    U0j(i) = p0j(i) * ((1/(D(j,i))).^(1/(m-1))/SUM0);
                end
                
                if sum(U0j.^m .* D(j,:)) + lambda * sum(p0j) > sum(Uopt.^m .* D(j,:)) + lambda * nnz(Uopt)
					checkrj = 0;
                else 
                    Uopt = U0j;
                    if sum(p0j) == 1 
                        checkrj = 0;
                    end
                end
            end

            U(j,:) = Uopt;
        end
        J(iter) = sum(sum(U.^m .* D)) + lambda * nnz(U>0);

		% Convergence check
		fprintf('Iteration %d: J = %.4f\n', iter, J(iter));
        % fprintf('Iteration %d: J = %.4f, lambda = %.4f\n', iter, J(iter), lambda);
		if iter > 1 && abs(J(iter) - J(iter-1)) < conv
            J = J(J>0);

            % Compute Xie & Beni index (Minimise)
            XieBeni = XB(X, U, P, m);
            fprintf('Xie and Beni Index: %.4f\n', XieBeni);

			break;
        end
    end

