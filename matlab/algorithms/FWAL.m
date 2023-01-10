function [X, constraintResidual] = FWAL(Q,A,ATranspose,b,beta0,maxit,verbose,dWaveMode,isSave,isVisualize)
%HCGM2 Summary of this function goes here
%   Detailed explanation goes here

n = size(Q,1);
if size(Q,2) ~= n, error('Q must be square.'); end
if (~exist('verbose', 'var')) verbose = 1; end
if (~exist('dWaveMode', 'var')) dWaveMode = 0; end
if (~exist('isSave', 'var')) isSave = 0; end
if (~exist('isVisualize', 'var')) isVisualize = 0; end

isEarlyStopping = 0;

tic; %start timer

X = zeros(n,n);
y = zeros(size(b));
constraintResidual = -b;
beta = beta0;

if (dWaveMode>0)
    %global dwaveModule;
    dwaveModule=load_dwave_module();
end

if (isVisualize)
    clrs = lines(7);
    markers = {'o', '+', '*', 'x', 'square', 'diamond', 'pentagram', 'hexagram', 'v'};
    figure, beautify_plot; grid on;
    constraintResiduals = zeros(maxit,1);
    costs = zeros(maxit,1);
    costsRounded = zeros(maxit,1);
    etas = zeros(maxit,1);
    AHbs = zeros(maxit,1);
    HQs = zeros(maxit,1);
end

minCost = 9999;
prevCost = 9999;
for t = 1:maxit
    
    % parameters
    beta = beta0 * sqrt(t+1); 
%     beta = beta * 1.01;
    eta = 2 / (t+1);
    
    % gradient evaluation
    grad = Q + ATranspose(y + beta .* constraintResidual);
    
    % linear minimization oracle
    if (dWaveMode==1)
        x = solve_dwave_module(dwaveModule, grad);
    elseif (dWaveMode==2)
        x = solve_dwave_module(dwaveModule, grad, 1);
    else
        x = solve_exhaustive_binary(grad);
    end    

    % convex combination step
    H = x*x';
    X =  X + eta .* (H - X);   
    
    constraintResidual = constraintResidual + eta * (A(x) - b - constraintResidual);
    
    if (isEarlyStopping)
        [qx, energyMunkres,~,~] = project_relaxed(grad, X);
        qx = [1; qx];
        H = qx*qx';
        if (Q(:)'*H(:) < Q(:)'*X(:))
            X = H;
            constraintResidual = A(qx)-b;
        end
    end
    
    % dual update
    y = y + beta0 * constraintResidual;

    curCost = Q(:)'*X(:);
    if (verbose)
        fprintf("Constraint: %2e, Cost: %2e \n",[norm(constraintResidual), curCost]);
    end
    if 2^floor(log2(t)) == t
        fprintf("Iteration: %d, Constraint: %2e, Cost: %2e \n",[t, norm(constraintResidual), curCost]);
    end

    if (isEarlyStopping)
        if (norm(constraintResidual)==0 && prevCost == curCost)
            break;
        end
        prevCost = curCost;
    end

    if isSave 
        % if mod(t,10) == 0
            if ~exist('saveDir','var')
                if (dWaveMode>0)
                    saveDir = ['results/Dwave_',num2str(dWaveMode),'_',datestr(now,30),'/'];
                else
                    saveDir = ['results/',datestr(now,30),'/'];
                end
                mkdir(saveDir);
            end
            time = toc;
            save([saveDir,num2str(t)]);
        % end
    end

    if (isVisualize)

        %HQs(t) = H(:)' * Q(:);
        %AHbs(t) = norm(A(x) - b);
        %plot(1:t, HQs(1:t),'LineWidth',3,'Marker',markers{2}, 'Color', clrs(1,:));
        %hold on, plot(1:t, AHbs(1:t),'LineWidth',3,'Marker',markers{3}, 'Color', clrs(2,:));
%
        %legend({'H(:)^T Q(:)', 'A(x)-b'});
        %drawnow;
        %hold off;
        %continue;

        % post-processs
        XFW = X(2:end, 2:end);
        [u,~] = svd(XFW);
        u = u(:,1);
        if (sum(u>0)<length(u)/2)
            u = - u;
        end

        PsX = perms_q_to_cell(u, 3);
        [PsX,~] = round_munkres(PsX);
        qx = perms_cell_to_q(PsX);
        curCostMunkres = qx'*Q(2:end,2:end)*qx;

        constraintError = norm(constraintResidual);
        constraintResiduals(t) = norm(constraintError);
        costs(t) = curCost;
        costsRounded(t) = curCostMunkres;
        plot(1:t, costs(1:t),'LineWidth',3,'Marker',markers{2}, 'Color', clrs(1,:));
        hold on, plot(1:t, constraintResiduals(1:t),'LineWidth',3,'Marker',markers{3}, 'Color', clrs(2,:));
        hold on, plot(1:t, costsRounded(1:t),'LineWidth',3,'Marker',markers{end}, 'Color', clrs(3,:));
        legend({'QUBO objective', 'constraint objective', 'rounded objective'});
        drawnow;
        hold off;
    end

%     if (curCost+norm(constraintResidual)<minCost)
%         minX = X; 
%         minConstraintResidual = constraintResidual;
%         minCost = curCost+norm(constraintResidual);
% %         if(t>8)
% %             break;
% %         end
%     end
% 
%     if 2^floor(log2(t)) == t
%         [uu,zz] = eig(X);
%         zz = diag(zz);
%         [~,ind] = max(abs(zz));
%         
%         AAA = perms_q_to_cell(uu(:,ind),2);
%         AAA{:}
%     end

end

clear dwaveModule;

end

function out = LMO(in)
    [~, out] = solve_exhaustive(grad, nCams, N);
end
