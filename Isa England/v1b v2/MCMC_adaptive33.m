function [xsto, outsto, history, accept_rate] = MCMC_adaptive33(F, x0, n, sigma, fixinds, blockind, cov0, displ)

d = length(x0); b = 0.05; sd = sigma*2.4^2/d; sd = 1;
if ~isempty(fixinds)
    inds = fixinds(1,:); vals = fixinds(2,:);
else
    inds = []; vals = [];
end

% Checks on the initial covariance matrix
if isempty(cov0)
    cov0 = eye(d); cov0(inds,:) = 0; cov0(:,inds) = 0;
end
cov0(1:blockind,blockind+1:end) = 0; cov0(blockind+1:end,1:blockind) = 0;

% Initiate the output matrices
xsto = zeros(d,n); outsto = zeros(1,n);
history = zeros(d+1,n);                                                    % Rows: 1:d Proposed values 4. Accept or reject

xsto(:,1) = x0(:); xbar = xsto;
FX = F(x0); outsto(1) = FX;
acc = 0;

if displ; figure; end

% --- Start the MCMC loop -------------------------------------------------
for t = 2:n
    try
    X = xsto(:,t-1);
    
    % --- Make a proposal from the distribution
    Y0 = mvnrnd(X,0.1^2*cov0*sigma/d);
    if t < 101
        Y = max(Y0,0); Y(inds) = vals;
    else
        ind0 = t-100; ind1 = t-1;
        covmat = cov(xsto(:,ind0:ind1)');
        covmat(inds,:) = 0; covmat(:,inds) = 0;
        covmat(1:blockind,(blockind+1:end)) = 0; covmat((blockind+1):end,1:blockind) = 0;
        % to make sure the matrix is positive semi def
        epsilon = 1e-6;  % the small regularization constant
        covmat = covmat + epsilon * eye(d);
        
        % Precaution to make sure values don't get negative
        Y = max((1-b)*mvnrnd(X,sd*covmat) + b*Y0,0);
        Y(inds) = vals;
    end
    history(1:d,t) = Y;
    
    % --- Decide whether to accept or not
    FY = F(Y);
    if (rand < exp(FY-FX)) && (abs(FY) < Inf)
        % Accept
        xsel = Y(:);
        FX = FY;
        acc = acc+1;
        history(end,t) = 1;
    else
        % Reject
        xsel = xsto(:,t-1);
    end
    xsto(:,t) = xsel;
    outsto(t) = FX;
    xbar(:,t) = (xbar(:,t-1)*(t-1) + xsel)/t;
    
    % Display options
    if displ && (mod(t,round(n/25))==0); fprintf('%0.5g ', t/n*25); end
    if displ && (mod(t,200)==0)
        plot(xsto(:,1:t-1)'); xlim([0 n]); drawnow;
    end
    catch
        keyboard;
    end
end

accept_rate = acc/n;
xsto = xsto';
history = history';