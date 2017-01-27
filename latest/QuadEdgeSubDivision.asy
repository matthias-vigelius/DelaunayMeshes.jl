import graph;

size(20cm,20cm);
draw((0,0)--(2,0), Arrow);
draw((2,0)--(1,2), Arrow);
draw((1,2)--(0,0), Arrow);
filldraw(Circle((0,0), 0.01));
filldraw(Circle((2,0), 0.01));
filldraw(Circle((1,2), 0.01));

label("$1$", (0,0), SW);
label("$2$", (2,0), SE);
label("$3$", (1,2), N);

draw(Circle((1,1), 0.01));
draw(Circle((1,-1), 0.01));

draw((1, -1)--(1,1), dashed, Arrow);
draw((1,1)..(2,1.5)..(3,0)..(1,-1), dashed, Arrow);
draw((1,1)..(0,1.5)..(-1,0)..(1,-1), dashed, Arrow);
