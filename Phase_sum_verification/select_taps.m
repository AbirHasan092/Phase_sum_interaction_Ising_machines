function sel = select_taps(offsets, M, mode, tol)
%SELECT_TAPS  Choose M physical taps near uniform phase positions and build
%   the complementary coupling map.
%
%   offsets : L x 1 measured tap offsets [rad], q_j(phi) ~ q_1(phi + offsets(j))
%   M       : number of coupling paths (2 <= M <= L)
%   mode    : 'exact'   - every target 2*pi*(r-1)/M must be matched within tol
%             'nearest' - closest distinct taps (greedy unique assignment)
%   tol     : tolerance [rad] for exact matching
%
%   Selected taps are ordered by target phase, r = 1..M (tap at target 0 first).
%   Complementary pairing  v_{1,r} <-> v_{2,c(r)},  c = [1, M, M-1, ..., 2]:
%       path r: receiving tap p(r) (oscillator 1), transmitting tap s(r) = p(c(r))
%   Since target(c(r)) = -target(r) mod 2*pi, the ideal arrangement is +delta_r / -delta_r.

offsets = offsets(:);
L = numel(offsets);
if M < 2 || M > L, error('select_taps: need 2 <= M <= L (M = %d, L = %d).', M, L); end
target = 2*pi*(0:M-1).'/M;
D = abs(angle(exp(1i*(offsets.' - target))));      % M x L circular distances

% exact uniform subset available?
[dmin, jmin] = min(D, [], 2);
sel.exact_available = all(dmin < tol) && numel(unique(jmin)) == M;

switch lower(mode)
    case 'exact'
        sel.available = sel.exact_available;
        p = jmin;
    case 'nearest'
        sel.available = true;
        p = zeros(M, 1);
        Dw = D;
        for it = 1:M                          % greedy: best remaining (target, tap) pair
            [~, lin] = min(Dw(:));
            [r, j] = ind2sub(size(Dw), lin);
            p(r) = j;
            Dw(r, :) = inf;  Dw(:, j) = inf;
        end
    otherwise
        error('select_taps: mode must be ''exact'' or ''nearest''.');
end
if numel(unique(p)) ~= M, error('select_taps: selected taps are not distinct.'); end

c = [1, M:-1:2];
sel.M = M;
sel.mode = lower(mode);
sel.p = p(:);
sel.s = p(c);
sel.s = sel.s(:);
sel.target = target;
sel.actual = offsets(p);
sel.dev = angle(exp(1i*(sel.actual - target)));                  % offset error per tap
sel.max_dev = max(abs(sel.dev));
sel.comp_err = angle(exp(1i*(offsets(sel.p) + offsets(sel.s)))); % ideal: delta_p + delta_s = 0
sel.max_comp_err = max(abs(sel.comp_err));
end
