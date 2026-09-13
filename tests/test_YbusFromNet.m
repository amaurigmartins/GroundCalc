function test_YbusFromNet()
% Run from the repository root: addpath('tests'); test_YbusFromNet
% No GUI, external data or additional toolboxes are required.
old_path = path;
cleanup = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

% A single shunt must not introduce a fictitious second branch.
check_close(YbusFromNet([1 1 0 4], 0, 1), 1/4);
check_close(YbusFromNet([1 0 1 4], [], 1), 1/4);
fprintf('PASS: single shunt, ground orientation and empty mutual list\n');

% Two parallel paths to ground add admittances, not impedances.
check_close(YbusFromNet([1 1 0 4; 2 0 1 4], 0, 1), 1/2);
fprintf('PASS: parallel shunts\n');

% A 2-ohm series branch feeding a 4-ohm shunt has input impedance 6 ohms.
net = [1 1 2 2; 2 2 0 4];
Y = YbusFromNet(net, 0, 2);
check_close(Y, [1/2 -1/2; -1/2 3/4]);
V = Y\[10; 0];
check_close(V, [60; 40]);
check_close(V(2)/4, 10); % Injected current equals current returned to earth.
check_close(YbusFromNet([1 2 1 2; 2 0 2 4], 0, 2), Y);
fprintf('PASS: series/shunt voltages, current balance and reversed branches\n');

% Complex impedances must retain phase and use a symmetric (not Hermitian)
% admittance matrix, as in the original real incidence-matrix formulation.
z1 = 2 + 3i;
z2 = 4 - 1i;
Y = YbusFromNet([1 1 2 z1; 2 2 0 z2], 0, 2);
check_close(Y, [1/z1 -1/z1; -1/z1 1/z1+1/z2]);
V = Y\[1; 0];
check_close(V, [z1+z2; z2]);
fprintf('PASS: complex impedances and voltages\n');

% Compare a larger uncoupled circuit against the primitive-matrix method.
n = 30;
net = [(1:n-1)' (1:n-1)' (2:n)' (2:n)'+1i];
net = [net; n 1 0 10; n+1 n 0 20];
A = zeros(n+1, n);
A(1:n-1,:) = diff(eye(n));
A(n,1) = 1;
A(n+1,n) = 1;
expected = A.' * diag(1./net(:,4)) * A;
check_close(YbusFromNet(net, 0, n), expected);
fprintf('PASS: primitive-matrix equivalence for a 30-node circuit\n');

% The coupled path must retain its existing behavior.
net = [1 1 0 2; 2 2 0 3];
check_close(YbusFromNet(net, [1 2 0.5], 2), [3 -0.5; -0.5 2]/5.75);
fprintf('PASS: nonzero mutual impedance\n');
fprintf('All YbusFromNet regression checks passed.\n');
end

function check_close(actual, expected)
assert(isequal(size(actual), size(expected)), 'Unexpected matrix size.');
assert(all(isfinite(actual(:))), 'Non-finite result.');
assert(norm(actual-expected, 'fro') <= 1e-11*max(1,norm(expected,'fro')), ...
    'Result differs from the reference circuit.');
end
