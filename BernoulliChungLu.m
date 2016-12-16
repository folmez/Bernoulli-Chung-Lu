function A = BernoulliChungLu(varargin)
% A = BERNOULLICHUNGLU(k) generates a random graph with the given degree
% sequence k. The algorithm is based on the Chung-Lu graph model.
% 
% The algorithm used here follows Section 1.2 of MLA Winlaw, M., H.
% DeSterck, and G. Sanders. An In-Depth Analysis of the Chung-Lu Model. No.
% LLNL--TR-678729. Lawrence Livermore National Lab.(LLNL), Livermore, CA
% (United States), 2015.
%
% Version 1.1 (Dec 2016)
%   - Excluding self-edges is added
%   - Input style changed
%
% Version 1.0 (Dec 2016)
%   - Initial version

% Input parameters
k = varargin{1};
nr_reps = 1;
exclude_self_edges = 0;

i = 2;
while i <= length(varargin),
    switch varargin{i},
        case 'nr_reps',             nr_reps = varargin{i+1};
        case 'exclude_self_edges',  exclude_self_edges = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

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
    if exclude_self_edges==1
        A = triu(U,1) < triu(P,1);
    else
        A = triu(U) < triu(P);
    end
    A = A+A';    
    
    % Calculate simulated degree sequence
    k_hat(:, i) = sum(A)';
end
k_hat_mean = mean(k_hat, 2);

% Plot degrees
figure,
plot(1:n, k, 'bo', 'MarkerSize', 10);
hold on;
if exclude_self_edges==1
    plot(1:n, k-k.^2/(2*m), 'ko', 'MarkerSize', 10);
end
plot(1:n, k_hat_mean, 'rx', 'MarkerSize', 10);
xlabel('Vertex');
ylabel('Degree');
if nr_reps==1 && exclude_self_edges==0
    h_leg = legend('Input degree', 'Simulation degree');
elseif nr_reps>1 && exclude_self_edges==0
    h_leg = legend('Input degree', ...
        ['Average over ' num2str(nr_reps) ' simulations']);
elseif nr_reps==1 && exclude_self_edges==1
    h_leg = legend('Input degree', 'Expected degree', 'Simulation degree');
elseif nr_reps>1 && exclude_self_edges==1
    h_leg = legend('Input degree', 'Expected degree', ...
        ['Average over ' num2str(nr_reps) ' simulations']);   
end
set(h_leg, 'Location', 'NorthWest');
xlim([0 n+1]);
ylim([min([k;k_hat_mean])-1 max([k;k_hat_mean])+1]);

end