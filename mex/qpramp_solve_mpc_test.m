timesteps = 100;

folder = "../../examples/example3";
a = csvread(folder + "/a.csv");
b = csvread(folder + "/b.csv");
x0 = csvread(folder + "/x0.csv");
invh = csvread(folder + "/invh.csv");
w = csvread(folder + "/w.csv");
g = csvread(folder + "/g.csv");
s = csvread(folder + "/s.csv");
f = csvread(folder + "/f.csv");

neg_g_invh_gt = -g*invh*g';
neg_s = -s;
neg_w = -w;
neg_invh_f = -invh*f;
neg_g_invh = -g*invh;

m = size(b, 2);
p = size(neg_invh_f, 1);

% Transposing in the loop gives big slowdown
neg_g_invh_gt_t = neg_g_invh_gt';
neg_s_t = neg_s';
neg_w_t = neg_w';
neg_invh_f_short_t = neg_invh_f(1:m,:)';
neg_g_invh_short_t = neg_g_invh(:,1:m)';

x = x0(1,:)';

tic;
for i = 1:timesteps
    u = qpramp_solve_mpc(p, neg_g_invh_gt_t, neg_s_t, neg_w_t, neg_invh_f_short_t, neg_g_invh_short_t, x');
    x = a*x + b*u';
end
toc;
disp(x');
