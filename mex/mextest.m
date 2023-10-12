m = 4;
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
neg_g_invh = -g*invh;

neg_invh_f = -invh*f;

x = x0(1,:)';

tic;
for i = 1:timesteps
    z = qp_ramp_solve(neg_g_invh_gt', neg_s', neg_w', neg_g_invh', x');
    u = z(1:m)' + neg_invh_f(1:m,:)*x;
    x = a*x + b*u;
end
toc;
disp(x');