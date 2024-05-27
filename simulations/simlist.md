# Run descriptions

1. White noise, c(1, 1, 1), interval 0, sight 1

2. White noise, c(1, 1, 1), interval 1, sight 1

3. (Same as 1) but par0 = 4, 0, -0.01

4. (Same as 1) but par_start = par0

== Change prep function ===

5. White noise, par_start = par0, interval 0, sight 1.

6. (Same as 5 except) par_start = c(1, 1, 1)

7. (Same as 6 except) interval = 1.

8. ^, interval = 2.

9. interval = 0, sight = 2.

10. interval = 0, sight = 3.

11. interval = 0, sight = 1, log_likelihood uses max value instead of taking value at true interval

12. replicating 11 with true interval just to check discrepancy (paths reused)

13. interval = 0, sight = 1, true interval, autocorrelated landscape.
