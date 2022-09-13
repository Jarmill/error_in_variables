# error_in_variables
Data-driven stabilization and robust control of linear systems corrupted by error in variables

"Data-Driven Superstabilizing Control of Error-in-Variables Discrete-Time Linear Systems"


## Experiments 
Algorithm 1
test_Full_SS.m      			full method for superstability, only works for 1st order system due to high complexity

Algorithm 2
test_Dual_SS.m   		           dual method for superstability, greatly reduce the complexity, works for 2nd order system


"Data-Driven Stabilizing and Robust Control of Discrete-Time Linear Systems with Error in Variables"

section 7
test_Dual_SS_all_noise.m            dual method considering all type of noises (dx, du, w), this covers test_Dual_SS.m

Algorithm 5,6,7
test_Dual_QS.m                            dual method for quadratic stability, worst case H2 works while Hinf currently has bug

Algorithm 5,6,7 + section 7
test_Dual_QS_all_noise.m            dual method considering all type of noises (dx, du, w), this covers test_Dual_QS.m


complexity.m                                  computes the number of constraints and variables
Dual_SS_manual.m                      manually defines the multiplication in the coefficients, it is much faster than Dual_SS.m while
					currently only supports SS


Contact Tianyu (dai.ti@northeastern.edu) for any questions regarding the code


