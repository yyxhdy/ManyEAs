package jmetal.metaheuristics.nsgaIII;

import java.util.HashMap;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.thetadea.ThetaDEA;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.ProblemFactory;
import jmetal.util.JMException;

public class NSGAIII_main {
	public static void main(String args[]) throws JMException, ClassNotFoundException{
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection; //Selection operator
		
		HashMap parameters; // Operator parameters
		
		Object[] params = { "Real", 9, 5};
		problem = (new ProblemFactory()).getProblem("DTLZ1", params);
		
		algorithm = new NSGAIII(problem);
		

		algorithm.setInputParameter("normalize", true);
		
		
		
		algorithm.setInputParameter("div1", 6);
		algorithm.setInputParameter("div2", 0);
		
		
		algorithm.setInputParameter("maxGenerations", 600);
		
		// Mutation and Crossover for Real codification
		parameters = new HashMap();
		parameters.put("probability", 1.0);
		parameters.put("distributionIndex", 30.0);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
				parameters);

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation",
				parameters);

		parameters = null;
		selection = SelectionFactory.getSelectionOperator("RandomSelection",
				parameters);

		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);
		
		SolutionSet population = algorithm.execute();
		
		population.printObjectivesToFile("FUN");
	}
}
