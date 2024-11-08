function [U, P, J, iter] = FKM_L0_P(X, C, m, lambda, conv, Max_iter, stand)

%% Standardization
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
Uperm=dec2bin(0:2^C-1) - '0';
Uperm(1,:)=[];
[K, ~] = size(Uperm);

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
    epsilon = 1e-10;
    for j = 1:N
	    for i = 1:C
            D(j,i) = max(sum((X(j,:)-P(i,:)).^2), epsilon);
    	end
    end
    % Update of the membership degrees
    SUM = sum((1./D).^(1/(m-1)),2);
	for j = 1:N
		fUjopt = Inf;
		for k=1:K
			if sum(Uperm(k,:))==1
				Utent=Uperm(k,:);
			else
				SUM0 = sum(Uperm(k,:).*(1./D(j,:)).^(1/(m-1)),2);
				for i = 1:C 
					Utent(i)=Uperm(k,i)*((1/(D(j,i))).^(1/(m-1))/SUM0);
				end
			end
			if sum(Utent.^m .* D(j, :)) + lambda * nnz(Utent) < fUjopt
				fUjopt = sum(Utent.^m .* D(j, :)) + lambda * nnz(Utent);
				Uopt = Utent;
			end
		end
		U(j,:) = Uopt;
	end
    % Convergence check
    J(iter) = sum(sum(U.^m .* D)) + lambda * nnz(U>0);
    if iter > 1 && abs(J(iter) - J(iter-1)) < conv
        J = J(J>0);
	    break;
    end
end
