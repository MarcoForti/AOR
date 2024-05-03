

function [U, P, J, iter] = FKM_L0_P(X, C, m, lambda, conv, Max_iter, stand)


%% standardization
if stand == 1
    Jm = eye(n) - (1/n)*ones(n);
	X = Jm*X/diag(std(X,1));
end

%% Initialization
[N, T] = size(X); 
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
		    % D(j,i) = sum((X(j,:)-P(i,:)).^2);
            D(j,i) = max(sum((X(j,:)-P(i,:)).^2), epsilon);
    	end
    end
    % Update of the membership degrees
    SUM = sum((1./D).^(1/(m-1)),2);

% % % Cycle for all the permutations % % % 

  for j = 1:N
        for i = 1:C 
		    U(j,i) = (1/(D(j,i))).^(1/(m-1))/SUM(j);
        end

    checkrj = 1;    
    Uopt = U(j,:);  
    U0j = U(j,:);   
    p0j = (U0j > 0);

    while checkrj == 1
        Uperm = zeros(1,C);
        V1 = ones(1, C);
        for ip = 1:(C-1)
            V1(ip) = 0;
            VP = perms(V1);
            Uperm = U0j.*VP;
            [righe_uniche, indici, ~] = unique(Uperm, 'rows', 'stable');  
            Uperm = Uperm(sort(indici), :);  
            [K, ~] = size(Uperm);
            for k = 1:K
                U0jP = Uperm(k,:);
                p0jP = (U0jP > 0);

                SUM0 = sum(p0jP.*(1./D(j,:)).^(1/(m-1)),2);
                for i = 1:C 
                    U0jP(i) = p0jP(i) * ((1/(D(j,i))).^(1/(m-1))/SUM0);
                end  % Fine ciclo C


                if sum(U0jP.^m .* D(j, :)) + lambda * nnz(p0jP) > sum(Uopt.^m .* D(j, :)) + lambda * nnz(Uopt)
                    checkrj = 0; 
                else 
                    Uopt = U0jP;
                    if sum(p0jP) == 1
                        checkrj = 0;
                    end
                end  

        
            end  
    
        end  
        
    end  
    U(j,:) = Uopt;

  end  

    % % % End while cycle % % % 
    J(iter) = sum(sum(U.^m .* D)) + lambda * nnz(U>0);
    % Convergence check
    fprintf('Iteration %d: J = %.4f\n', iter, J(iter));
    if iter > 1 && abs(J(iter) - J(iter-1)) < conv
        J = J(J>0);
	    break;  
    end
end  




