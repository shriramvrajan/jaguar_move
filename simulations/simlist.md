(Change in prep function - cellFromXY and cellFromRowCol do NOT do the same thing!)

5. White noise, par_start = par0, interval 0, sight 1.

6. (Same as 5 except) par_start = c(1, 1, 1)

7. (Same as 6 except) interval = 1.

8. ^, interval = 2.

9. interval = 0, sight = 2.

10. interval = 0, sight = 3.

11. interval = 0, sight = 1, log_likelihood uses max value instead of taking value at true interval

12. replicating 11 with true interval just to check discrepancy (paths reused)

13. interval = 0, sight = 1, true interval, autocorrelated landscape.

14. interval = 0, sight = 1, true interval, autocorrelated landscape, non-monotonic env function

15. (Same as 14 except) white noise landscape

16all. Same as 6 except fitting all together instead of individually.

17all. Same as 14 except fitting all together. 

18all. Same as 17 except starting at par0.

19. Same as 17 and 18 except starting at 17 par_out.

20. Same as 17 except sim_steps not fixed at true value.

21. Same as 17 except not monotonic; and env normalized.

22. Not monotonic; env normalized; fit individually.
