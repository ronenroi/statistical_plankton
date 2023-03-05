close all
%% clerr
errors = commonline_errors(clerr(:,3),indeces,size(projections,3));
[B,I] = sort(errors);
figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
figure;plot(B)
%% Pij
errors = commonline_errors(Pij,indeces,size(projections,3));
[B,I] = sort(errors);
figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
figure;plot(B)
%% W based
COV = W +eye(size(W,1));
    
    PC = pcacov(COV);
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1);
    
    [B,I] = sort(valid_eq);
figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
figure;plot(B)
%% W norm based
% Prepare weights for the blocks of the matrix S.
D = mean(W,2); % We use mean and not sum so that after normalizing W the 
               % sum of its rows will be N, as required.
nulls = abs(D)<1.0e-13;
Dhalf = D;
Dhalf(~nulls) = Dhalf(~nulls).^-0.5;
Dhalf(nulls) = 0;
Dhalf = diag(Dhalf);

% Assert Rows Normalization
W_normalized = (Dhalf.^2)*W;
COV = W_normalized +eye(size(W_normalized,1));
    
    PC = pcacov(COV);
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1);
    
    [B,I] = sort(valid_eq);
figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
figure;plot(B)
%% A
COV = A + A' ;
    
    PC = pcacov(COV);
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1);
    [B,I] = sort(valid_eq);
figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
figure;plot(B)
%% correlation refined based
COV = corrstack_refined + corrstack_refined' + eye(size(corrstack_refined,1));
    
    PC = pcacov(COV);
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1);
    [B,I] = sort(valid_eq);
figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
figure;plot(B)
%% correlation based
COV = corrstack + corrstack' + eye(size(corrstack,1));
    
    PC = pcacov(COV);
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1);
    [B,I] = sort(valid_eq);
figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
figure;plot(B)