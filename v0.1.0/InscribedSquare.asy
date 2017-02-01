import graph;

size(20cm,20cm);

pair v1 = (1,1);
pair v2 = (2,1);
pair v3 = (1.5,2);
draw(v1 -- v2 -- v3 -- v1);

real dx = 0.5;
pair sc = (1.5, 1.25);
pair s1 = sc + 0.5 * dx * (S + W);
pair s2 = sc + 0.5 * dx * (S + E);
pair s3 = sc + 0.5 * dx * (N + E);
pair s4 = sc + 0.5 * dx * (N + W);

draw(s1 -- s2 -- s3 -- s4 -- s1, dashed);

pair ss1 = sc + 0.05 * dx * (S + W);
pair ss2 = sc + 0.05 * dx * (S + E);
pair ss3 = sc + 0.05 * dx * (N + E);
pair ss4 = sc + 0.05 * dx * (N + W);

draw(ss1 -- ss2 -- ss3 -- ss4 -- ss1, blue + linewidth(1mm));
