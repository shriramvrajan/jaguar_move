# Simulations

1 corresponds to sim1 folder, etc.

1. 2 parameters, simple exponential functional form, with parameter landscapes

2. 4 parameters, 2nd order polynomial function, poor resolution because of bad par0 values.

3. 4 parameters, 2nd order polynomial, better resolution. Plotted generating function vs fitted functions and there is systematic bias in the estimation. The stay-move parameter values are accurate, not to the generating value (~0.88) but to the actually generated frequencies (~0.75).

4. 3 parameters, 2nd order polynomial, removed movement kernels. Fits get worse. Still have to check the tradSSF for both the previous sim and this one.

5. 4 parameters, 2op, with mks. But I **fixed the move/stay bug**. Movement parameter of 2, i.e. Pm of ~0.89.

6. 4p, 2op, movement parameter of -2, i.e. Pm of ~0.11.

7. 4p, 2op, Pm of ~0.11, 100 runs.

8. 4p, 2op, Pm of ~0.88, 50 runs.

9. 4p, 2op, Pm of ~0.88, white noise landscape.

10. 4p, 2op, truer white noise landscape.

11. Uniform landscape - all raster cells = 1.

12. Runif landscape - unif dist from 0 to 8.

13. Uniform landscape, more trials.

14. New method: not normalizing twice. White noise landscape.

15. Replicate of 14.
