package jmetal.util.ranking;

import java.util.ArrayList;

import java.util.List;




import jmetal.core.Solution;
import jmetal.core.SolutionSet;

import jmetal.util.comparators.FitnessComparator;




public class ThetaRanking implements Ranking {

	private SolutionSet solutionSet_;

	private List<SolutionSet> ranking_;

	private SolutionSet[] refSets_;
	

	double[][] lambda_;
	double[] zideal_;

	double theta_;

	int obj_;
	
	final double inf = 1E6;
	
	boolean normalize_;

	public ThetaRanking(SolutionSet solutionSet,  double[][] lambda, 
			double[] zideal, 
			double theta, boolean normalize) {
		this.solutionSet_ = solutionSet;
		this.lambda_ = lambda;
	
		this.theta_ = theta;
		this.obj_ = solutionSet.get(0).numberOfObjectives();
	
		this.zideal_ = zideal;
		
		this.normalize_ = normalize;

		ranking_ = new ArrayList<SolutionSet>();

		refSets_ = new SolutionSet[lambda_.length];
		for (int i = 0; i < refSets_.length; i++)
			refSets_[i] = new SolutionSet();
		
		associate();
		rank();
		
	}
	
	
	
	void associate() {
		
		double[] dists = null;
		for (int k = 0; k < solutionSet_.size(); k++) {

			Solution sol = solutionSet_.get(k);

			dists = getDistances(sol, lambda_[0]);
			double d2 = dists[1];
			double d1 = dists[0];
			int index = 0;

			for (int j = 1; j < lambda_.length; j++) {

				dists = getDistances(sol, lambda_[j]);

				if (dists[1] < d2) {
					d2 = dists[1];
					d1 = dists[0];
					index = j;
				}
			}

			
			setFitness(sol, index, d1, d2);
			refSets_[index].add(sol);
			sol.setClusterID(index);

		}
	}
	
	
	void setFitness(Solution sol, int index, double d1, double d2){
		if (this.normalize_) {
			if (!isObjAxis(index))
				sol.setFitness(d1 + theta_ * d2);
			else
				sol.setFitness(d1 + inf * d2);
		}
		else
			sol.setFitness(d1 + theta_ * d2);
	}
	
	double[] getDistances(Solution sol, double[] ref){
		if (this.normalize_)
			return getDistancesWithNormalize(sol, ref);
		else
			return getDistancesWithoutNormalize(sol, ref);
	}
	
	
	
	boolean isObjAxis(int index){
		for (int j = 0; j < obj_; j++){
			if (lambda_[index][j] != 0 && lambda_[index][j] != 1)
				return false;
		}
		return true;
	}
	
	
	
	double[] getDistancesWithoutNormalize(Solution sol, double[] ref) {
		double[] d = new double[2];

		double d1, d2, nl;

		d1 = d2 = nl = 0.0;

		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			d1 += (sol.getObjective(i) - zideal_[i]) * ref[i];
			nl += (ref[i] * ref[i]);
		}
		nl = Math.sqrt(nl);
		d1 = Math.abs(d1) / nl;

		d2 = 0;
		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			d2 += ((sol.getObjective(i) - zideal_[i]) - d1 * (ref[i] / nl))
					* ((sol.getObjective(i) - zideal_[i]) - d1 * (ref[i] / nl));
		}
		d2 = Math.sqrt(d2);

		d[0] = d1;
		d[1] = d2;

		return d;
	}

	double[] getDistancesWithNormalize(Solution sol, double[] ref) {
		double ip = 0;
		double refLenSQ = 0;

		double[] d = new double[2];
		for (int j = 0; j < obj_; j++) {

			ip += sol.getNormalizedObjective(j) * ref[j];
			refLenSQ += (ref[j] * ref[j]);
		}
		refLenSQ = Math.sqrt(refLenSQ);

		d[0] = Math.abs(ip) / refLenSQ;

		d[1] = 0;
		
		
		
		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			d[1] += (sol.getNormalizedObjective(i) - d[0] * (ref[i] / refLenSQ))
					* (sol.getNormalizedObjective(i) - d[0]
							* (ref[i] / refLenSQ));
		}
		d[1] = Math.sqrt(d[1]);
		
		return d;
	}


	
	

	
	
	
	void rank() {

		int maxLen = Integer.MIN_VALUE;
		for (int i = 0; i < refSets_.length; i++) {
			
			if (refSets_[i].size() > maxLen)
				maxLen = refSets_[i].size();
			refSets_[i].sort(new FitnessComparator());
		}


		for (int i = 0; i < maxLen; i++) {
			SolutionSet set = new SolutionSet();
			for (int j = 0; j < refSets_.length; j++) {
				if (refSets_[j].size() > i) {
					refSets_[j].get(i).setRank(i);
					set.add(refSets_[j].get(i));
				}
			}
			ranking_.add(set);
		}
	}


	
	
	public SolutionSet getSubfront(int rank) {
		return ranking_.get(rank);
	} // getSubFront

	/**
	 * Returns the total number of subFronts founds.
	 */
	public int getNumberOfSubfronts() {
		return ranking_.size();
	} // getNumberOfSubfronts
	
}
