

function[U,P,J] = FKM(X,m,C,conv)

[N,T] = size(X); 

%% Initialization

Y = rand(N,C);  

for j = 1:N
    U(j,:) = Y(j,:)/sum(Y(j,:));  
end
for i = 1:C
    for s=1:T
        P(i,s) = ((U(:,i).^m)'*X(:,s))/sum((U(:,i).^m)); 
    end
end
for j = 1:N
    for i = 1:C
	    D(j,i) = sum((X(j,:) - P(i,:)).^2);  
    end
end


%% Optimization

iter = 0;
Uold = U + 1;
while sum(sum((Uold-U).^2)) > conv  
    Uold = U;
    iter = iter+1;
    for i = 1:C
        for s = 1:T
            P(i,s) = ((U(:,i).^m)'*X(:,s))/sum((U(:,i).^m)); 
        end
    end
    for j = 1:N
        for i = 1:C
            D(j,i) = sum((X(j,:)-P(i,:)).^2);
        end
    end
    for j = 1:N
        SUM(j) = 0;
        for i = 1:C
            SUM(j) = SUM(j) + (1/(D(j,i))).^(1/(m-1));
        end
    end
    for j = 1:N
        for i = 1:C
            U(j,i) = (1/(D(j,i))).^(1/(m-1))/SUM(j);  
        end
    end
end

J = sum(sum(U.^m .* D));

