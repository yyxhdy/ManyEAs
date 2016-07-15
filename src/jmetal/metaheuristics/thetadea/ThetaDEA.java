package jmetal.metaheuristics.thetadea;



import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import Jama.Matrix;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;

import jmetal.operators.mutation.Mutation;

import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Permutation;
import jmetal.util.PseudoRandom;

import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.FitnessComparator;

import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.ranking.ThetaRanking;

import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;


public class ThetaDEA extends Algorithm {
	
	private int populationSize_;   // population size
	

	private SolutionSet population_;   // current population
	SolutionSet offspringPopulation_;  // offspring population
	
	SolutionSet union_;    // the union of current population and offspring population
	
	
	
	int generations_;   // generations

	
	/* if only one layer is adopted, div2_=0 */
	int div1_;  // divisions in the boundary layer
	int div2_;  // divisions in the inside layer
	
	double theta_;     // parameter theta 
	
	Operator crossover_; // crossover
	
	
	Operator mutation_;   // mutation operator


	
	boolean normalize_;  // normalization or not
	

	double[][] lambda_; // reference points
	

	double[] zideal_;   // ideal point

	double[] znadir_;   // nadir point
	
	double[][] extremePoints_; // extreme points
	
	
	public ThetaDEA(Problem problem) {
		super(problem);
	} // ThetaDEA 
	

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxGenerations;  // maximum number of generations

		generations_ = 0;
		
		
		/* set parameters */
		maxGenerations = ((Integer) this.getInputParameter("maxGenerations"))
				.intValue();

		theta_ =  ((Double)this.getInputParameter("theta")).doubleValue();

		div1_ = ((Integer) this.getInputParameter("div1")).intValue();

		div2_ = ((Integer) this.getInputParameter("div2")).intValue();

		normalize_ = ((Boolean) this.getInputParameter("normalize"))
				.booleanValue();
		

		/* generate two-layer weight vectors */
		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,
				problem_.getNumberOfObjectives());
		lambda_ = vg.getVectors();

		
		
		/*the population size is the same with the number of weight vectors*/
		populationSize_ = vg.getVectors().length; 
		

		crossover_ = operators_.get("crossover"); // set the crossover operator
		mutation_ = operators_.get("mutation");  // set the mutation operator
		
	
		initPopulation();   // initialize the population;
		
		initIdealPoint();  // initialize the ideal point
		
		initNadirPoint();    // initialize the nadir point
		
		initExtremePoints(); // initialize the extreme points
		
		
		while (generations_ < maxGenerations) {
			
			createOffSpringPopulation();  // create the offspring population
			
			union_ = population_.union(offspringPopulation_);
	
			SolutionSet[] sets = getParetoFronts();
			
			SolutionSet firstFront = sets[0];   // the first non-dominated front
			SolutionSet stPopulation = sets[1]; // the population used in theta-non-dominated ranking
			
			updateIdealPoint(firstFront);  // update the ideal point
			
			if (normalize_){
				updateNadirPoint(firstFront);  // update the nadir point
				normalizePopulation(stPopulation);  // normalize the population using ideal point and nadir point
			}

			getNextPopulation(stPopulation);  // select the next population using theta-non-dominated ranking
			
			generations_++;
			
		}
		
		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);
		
	}
	
	
	public void initExtremePoints() {
		int obj = problem_.getNumberOfObjectives();
		extremePoints_ = new double[obj][obj];
		for (int i = 0; i < obj; i++){
			for (int j = 0; j < obj; j++){
				extremePoints_[i][j] = 1.0e+30;
			}
		}
		
	}
	
	void getNextPopulation(SolutionSet pop){
		Ranking ranking = new ThetaRanking(pop, lambda_, zideal_, 
				theta_, normalize_);
		
		int remain = populationSize_;
		int index = 0;
		SolutionSet front = null;
		population_.clear();

		// Obtain the next front
		front = ranking.getSubfront(index);

		
		while ((remain > 0) && (remain >= front.size())) {
			
			
			for (int k = 0; k < front.size(); k++) {
				population_.add(front.get(k));
			} // for

			// Decrement remain
			remain = remain - front.size();

			// Obtain the next front
			index++;
			if (remain > 0) {
				front = ranking.getSubfront(index);
			} // if
		} // while

		if (remain > 0) { // front contains individuals to insert
		
			int[] perm = new Permutation().intPermutation(front.size());
			for (int k = 0; k < remain; k++) {
				population_.add(front.get(perm[k]));
			} // for
			remain = 0;
			
		} // if
	}

	SolutionSet[] getParetoFronts() {
		
		SolutionSet[] sets = new SolutionSet[2];
		Ranking ranking = new NondominatedRanking(union_);

		int remain = populationSize_;
		int index = 0;
		SolutionSet front = null;
		SolutionSet mgPopulation = new SolutionSet();

		front = ranking.getSubfront(index);

		sets[0] = front;

		while ((remain > 0) && (remain >= front.size())) {

			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			} // for

			// Decrement remain
			remain = remain - front.size();

			// Obtain the next front
			index++;
			if (remain > 0) {
				front = ranking.getSubfront(index);
			} // if
		}
		if (remain > 0) { // front contains individuals to insert
			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			}
		}

		sets[1] = mgPopulation;

		return sets;
	}

	void initPopulation() throws JMException, ClassNotFoundException {
		
		population_= new SolutionSet(populationSize_);
		
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population_.add(newSolution);
		}
	} 

	
	void createOffSpringPopulation() throws JMException {
		offspringPopulation_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) 
			doCrossover(i);
	}
	
	
	void doCrossover(int i) throws JMException{
		int r;
		do {
			r = PseudoRandom.randInt(0, populationSize_ - 1);
		} while (r == i);
		
		Solution[] parents = new Solution[2];
		
		parents[0] = population_.get(i);
		parents[1] = population_.get(r);
		
		Solution[] offSpring = (Solution[]) crossover_
				.execute(parents);
		
		mutation_.execute(offSpring[0]);
		
		problem_.evaluate(offSpring[0]);
		
		
		offspringPopulation_.add(offSpring[0]);
	}
	
	
	void copyObjectiveValues(double[] array, Solution individual) {
		for (int i = 0; i < individual.numberOfObjectives(); i++) {
			array[i] = individual.getObjective(i);
		}
	}
	

	double asfFunction(Solution sol, int j) {
		double max = Double.MIN_VALUE;
		double epsilon = 1.0E-6;

		int obj = problem_.getNumberOfObjectives();

		for (int i = 0; i < obj; i++) {

			double val = Math.abs((sol.getObjective(i) - zideal_[i])
					/ (znadir_[i] - zideal_[i]));

			if (j != i)
				val = val / epsilon;

			if (val > max)
				max = val;
		}

		return max;
	}

	double asfFunction(double[] ref, int j) {
		double max = Double.MIN_VALUE;
		double epsilon = 1.0E-6;

		int obj = problem_.getNumberOfObjectives();

		for (int i = 0; i < obj; i++) {

			double val = Math.abs((ref[i] - zideal_[i])
					/ (znadir_[i] - zideal_[i]));
			

			if (j != i)
				val = val / epsilon;

			if (val > max)
				max = val;
		}

		return max;
	}
	
	
	
	void initIdealPoint() {
		int obj = problem_.getNumberOfObjectives();
		zideal_ = new double[obj];
		for (int j = 0; j < obj; j++) {
			zideal_[j] = Double.MAX_VALUE;

			for (int i = 0; i < population_.size(); i++) {
				if (population_.get(i).getObjective(j) < zideal_[j])
					zideal_[j] = population_.get(i).getObjective(j);
			}
		}
	}
	
	
	void updateIdealPoint(SolutionSet pop){
		for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
			for (int i = 0; i < pop.size(); i++) {
				if (pop.get(i).getObjective(j) < zideal_[j])
					zideal_[j] = pop.get(i).getObjective(j);
			}
		}
	}
	
	void initNadirPoint() {
		int obj = problem_.getNumberOfObjectives();
		znadir_ = new double[obj];
		for (int j = 0; j < obj; j++) {
			znadir_[j] = Double.MIN_VALUE;

			for (int i = 0; i < population_.size(); i++) {
				if (population_.get(i).getObjective(j) > znadir_[j])
					znadir_[j] = population_.get(i).getObjective(j);
			}
		}
	}
	
	
	

	
	void updateNadirPoint(SolutionSet pop){
		
		updateExtremePoints(pop);

		
		int obj = problem_.getNumberOfObjectives();
		double[][] temp = new double[obj][obj];

		for (int i = 0; i < obj; i++) {
			for (int j = 0; j < obj; j++) {
				double val = extremePoints_[i][j] - zideal_[j];
				temp[i][j] = val;
			}
		}

		Matrix EX = new Matrix(temp);

		boolean sucess = true;
		
		if (EX.rank() == EX.getRowDimension()) {
			double[] u = new double[obj];
			for (int j = 0; j < obj; j++)
				u[j] = 1;

			Matrix UM = new Matrix(u, obj);

			Matrix AL = EX.inverse().times(UM);

			int j = 0;
			for (j = 0; j < obj; j++) {

				double aj = 1.0 / AL.get(j, 0) + zideal_[j];
		

				if ((aj > zideal_[j]) && (!Double.isInfinite(aj)) && (!Double.isNaN(aj)))
					znadir_[j] = aj;
				else {
					sucess = false;
					break;
				}
			}
		} 
		else 
			sucess = false;
		
		
		if (!sucess){
			double[] zmax = computeMaxPoint(pop);
			for (int j = 0; j < obj; j++) {
				znadir_[j] = zmax[j];
			}
		}
	}
	
	
	
	
	public void updateExtremePoints(SolutionSet pop){
		for (int i = 0; i < pop.size(); i++)
			updateExtremePoints(pop.get(i));
	}
	
	
	public void updateExtremePoints(Solution individual){
		int obj = problem_.getNumberOfObjectives();
		for (int i = 0; i < obj; i++){
			double asf1 = asfFunction(individual, i);
			double asf2 = asfFunction(extremePoints_[i], i);
			
			if (asf1 < asf2){
				copyObjectiveValues(extremePoints_[i], individual);
			}
		}
	}
	
	
	double[] computeMaxPoint(SolutionSet pop){
		int obj = problem_.getNumberOfObjectives();
		double zmax[] = new double[obj];
		for (int j = 0; j < obj; j++) {
			zmax[j] = Double.MIN_VALUE;

			for (int i = 0; i < pop.size(); i++) {
				if (pop.get(i).getObjective(j) > zmax[j])
					zmax[j] = pop.get(i).getObjective(j);
			}
		}
		return zmax;
	}
	
	void normalizePopulation(SolutionSet pop) {

		int obj = problem_.getNumberOfObjectives();

		for (int i = 0; i < pop.size(); i++) {
			Solution sol = pop.get(i);

			for (int j = 0; j < obj; j++) {

				double val = (sol.getObjective(j) - zideal_[j])
						/ (znadir_[j] - zideal_[j]);

				sol.setNormalizedObjective(j, val);
			}
		}
	}
	

}
