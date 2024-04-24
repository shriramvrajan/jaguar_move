# List of simulations

1. Low movement. Good parameter predictions for median but not mean.

2. Low movement, but with sim_interval fixed.

3. Sim_interval fixed, no movement kernel, only environment. Median works but not mean.

4. Replicate of 3. 

5. Trying out 25-cell neighborhood instead of 9.

6. Trying out 49-cell neighborhood.

7. Brian's suggestion: 9cell, white noise, sim_interval = step_size. PAR_FIT = PAR_TRUE

8. sim_interval = 2, white noise

9. sim_interval = 2, step_size = 2, white noise

10. sim_interval = 1, step_size = 2, white noise

11. sim_interval = 1, step_size = 1, SAC landscape

12. sim_interval = 1, step_size = 3, SAC landscape

13. sim_interval = 2, step_size = 1, SAC landscape

14. sim_interval = 2, step_size = 3, SAC landscape

15. replicate 11 but starting from 0,0,0 instead of par0. Optim doesn't move away from 0.

16. Replicate 11 starting from 5,5,5.

17. replicate 11 starting from 1,1,1.

18. FIXED SIM INTERVAL BUG. Trying sim_interval = 1, step_size = 1.

19. sim_interval = 2, step_size = 1.
