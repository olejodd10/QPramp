input_folder = "../data/model3";
initial_condition = 0;
timesteps = 100;

a = csvread(input_folder + "/a.csv");
b = csvread(input_folder + "/b.csv");
x0 = csvread(input_folder + "/x0.csv");
invh = csvread(input_folder + "/invh.csv");
w = csvread(input_folder + "/w.csv");
g = csvread(input_folder + "/g.csv");
s = csvread(input_folder + "/s.csv");
f = csvread(input_folder + "/f.csv");

neg_g_invh_gt = -g*invh*g';
neg_s = -s;
neg_w = -w;
neg_g_invh = -g*invh;

m = size(b, 2);
neg_invh_f = -invh*f;

% Transposing in the loop gives big slowdown
neg_g_invh_gt_t = neg_g_invh_gt';
neg_s_t = neg_s';
neg_w_t = neg_w';
neg_g_invh_t = neg_g_invh';

x = x0(initial_condition+1,:)';

tic;
for i = 1:timesteps
    z = qpramp_solve(neg_g_invh_gt_t, neg_s_t, neg_w_t, neg_g_invh_t, x');
    u = z(1:m)' + neg_invh_f(1:m,:)*x;
    x = a*x + b*u;
end
t = toc;
fprintf("Total time for initial condition %d running %d timesteps: %.0f us\n", initial_condition, timesteps, t*1e6);
disp(x');
clear mex; % Free memory used by solver
