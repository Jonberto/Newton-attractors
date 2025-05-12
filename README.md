This project visualizes the behavior of Newton's method applied to complex polynomials of the form f(x) = x^d - 1, using parallel computation in C with the C11 thread library. It generates two PPM images: one showing the attractors (which root each point converges to), and one showing the convergence speed (number of iterations until convergence). Each pixel corresponds to a point in the complex plane, and all iterations are computed independently which allows us to parallelize the computation efficiently.

Ex to run it ./newton -t5 -l1000 7
-t number of threads, -l specifies the resolution of the ppm image and argument "7" the degree of the polynomial 
