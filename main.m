% MIS enumeration for the union confusion graph (q=3, Xi in {0,1,2})
%   s ~ t  iff  Xi(s) == Xi(t)  AND  fi(s) ~= fi(t)
% method of Tsukiyama
clear; clc;

% ---- Vertex set: states (x1,x2,x3) in {0,1,2}^3 ----
states = zeros(27,3);
labels = cell(27,1);
idx = 1;
for x1 = 0:2
    for x2 = 0:2
        for x3 = 0:2
            states(idx,:) = [x1,x2,x3];
            labels{idx} = sprintf('(%d,%d,%d)', x1,x2,x3);
            idx = idx + 1;
        end
    end
end
n = size(states,1);

% ---- Functions (mod 3) ----
f1 = @(x) mod(x(2) + x(3), 3); %   f1 = x2 + x3 (mod 3)
f2 = @(x) mod(x(1) + x(3), 3); %   f2 = x1 + x3 (mod 3)
f3 = @(x) mod(x(1) + 2*x(2), 3);    %   f3 = x1 + 2*x2 (mod 3)

% ---- Build union adjacency ----
A = false(n,n);
for i = 1:n
    xi = states(i,:);
    fi = [f1(xi), f2(xi), f3(xi)];
    for j = i+1:n
        xj = states(j,:);
        fj = [f1(xj), f2(xj), f3(xj)];
        e1 = (xi(1) == xj(1)) && (fi(1) ~= fj(1));  % receiver 1
        e2 = (xi(2) == xj(2)) && (fi(2) ~= fj(2));  % receiver 2
        e3 = (xi(3) == xj(3)) && (fi(3) ~= fj(3));  % receiver 3
        A(i,j) = e1 || e2 || e3;
        A(j,i) = A(i,j);
    end
end
A(1:n+1:end) = false;           %zero diagonal
assert(isequal(A,A.'), 'Adjacency must be symmetric');

% ---- Enumerate MIS as cliques of the complement ----
Gc = ~A; Gc(1:n+1:end) = false; % complement without loops
MIS = bronKerboschPivot(Gc);    % each cell is a vector of vertex indices

% ---- Output ----
fprintf('Total maximal independent sets: %d\n\n', numel(MIS));  % ==> 225
for k = 1:numel(MIS)
    S = sort(MIS{k});
    fprintf('MIS %d (size %d): ', k, numel(S));
    fprintf('%d ', S); fprintf('\n    States: ');
    for t = 1:numel(S), fprintf('%s ', labels{S(t)}); end
    fprintf('\n');
end
%reza..checkpoint
% (Optional) size histogram
sizes = cellfun(@numel, MIS);
tab = tabulate(sizes); 
disp('Size distribution [size, count, percent]:');
disp(tab);

% Helper: Bron–Kerbosch with pivot ====
function MIS = bronKerboschPivot(G)
    n = size(G,1);
    MIS = {};
    MIS = bk([], 1:n, [], G, MIS);
end

function MIS = bk(R, P, X, G, MIS)
    if isempty(P) && isempty(X)
        MIS{end+1} = R; return;
    end
    % Choose pivot u maximizing |P ∩ N(u)|
    UX = [P, X];
    if isempty(UX)
        Nu = false(1,size(G,1));
    else
        best = -inf; uBest = UX(1);
        for u = UX
            cand = sum(G(u,P));
            if cand > best, best = cand; uBest = u; end
        end
        Nu = G(uBest,:);
    end
    % Iterate over P \ N(u)
    candidates = setdiff(P, find(Nu));
    Pset = P; Xset = X;
    for v = candidates
        Nv = find(G(v,:));
        MIS = bk([R v], intersect(Pset, Nv), intersect(Xset, Nv), G, MIS);
        Pset = setdiff(Pset, v);   % move v from P to X
        Xset = union(Xset, v);
    end
end
%mrds..checkpoint
