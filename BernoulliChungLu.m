function A = BernoulliChungLu(k, nr_reps)
% The algorithm used here follows Section 1.2 of MLA Winlaw, M., H.
% DeSterck, and G. Sanders. An In-Depth Analysis of the Chung-Lu Model. No.
% LLNL--TR-678729. Lawrence Livermore National Lab.(LLNL), Livermore, CA
% (United States), 2015.
%
% Version 1.0: December 2016

% Check inputs
if mod(sum(k),2)~=0
    error('Sum of degrees must be even!');
end

% Calculate number of edges
m = sum(k)/2;

% Check inputs
if any(k.^2>2*m)
    error('Square of every degree must be <= the number of edges');
end

% Find number of vertices
n = length(k);

% Reshape inputted degree sequence to a column vector
k = reshape(k, n, 1);

% Calculate edge probabilities
P = (k*ones(1,n)) .* (ones(n,1)*k');
P = P/(2*m);
P(1:n+1:n^2) = P(1:n+1:n^2)/2; % self edges are counted twice!

k_hat = zeros(n, nr_reps);
for i=1:nr_reps
    % Generate uniformly distributed nxn matrix
    U = rand(n);
    
    % Calculate adjacency matrix
    A = triu(U)<triu(P);
    A = A+A';
    
    % Calculate simulated degree sequence
    k_hat(:, i) = sum(A)';
end
k_hat_mean = mean(k_hat, 2);

% Plot degrees
figure,
plot(1:n, k, 'bo', 'MarkerSize', 10);
hold on;
plot(1:n, k_hat_mean, 'rx', 'MarkerSize', 10);
xlabel('Vertex');
ylabel('Degree');
if nr_reps==1
    legend('Input', 'Simulation', 'Location', 'North');
else
    legend('Input', ['Average over ' num2str(nr_reps) ' simulations'], ...
        'Location', 'North');
end
xlim([0 n+1]);
ylim([min([k;k_hat_mean])-1 max([k;k_hat_mean])+1]);

end