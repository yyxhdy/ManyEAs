# ManyEAs

This program mainly implements the Theta-DEA algorithm which was described in our paper:

Yuan Yuan, Hua Xu, Bo Wang, and Xin Yao. A new dominance relation-based evolutionary algorithm 
for many-objective optimization. 
IEEE Transactions on Evolutionary Computation, 2016, 20(1): 16-37.


It also includes an ''unofficial'' implementation of NSGA-III which was described in following paper:

Deb Kalyanmoy and Himanshu Jain. An evolutionary many-objective optimization algorithm using reference-point-based 
nondominated sorting approach, part I: Solving problems with box constraints. 
IEEE Transactions on Evolutionary Computation, 2014, 18(4): 577-601.


**************************************************************************************************************************************

The code was written in java and based on the jMetal framework. 
Moreover, it used an external archive, i.e., Jama-1.0.2.jar, which 
can be downloaded from http://math.nist.gov/javanumerics/jama/, 
so make sure you have linked it to the build path before you run the program. 
(In the current version, we used JAMA to compute the inverse matrix.)


The algorithmic procedures of Theta-DEA and NSGA-III can be seen in the files 
ThetaDEA.java and NSGAIII.java. 
We also illustrated how to run 
the two algorithms on a specific benchmark problem in the files ThetaDEA_main.java and NSGAIII_main.java, respectively.

Note that, this code can be used only for non-commercial purposes. 
We'd appreciate your acknowledgement if you use the code. 

For any problem concerning the code, please feel free to contact Dr. Yuan Yuan (yyxhdy@gmail.com).
