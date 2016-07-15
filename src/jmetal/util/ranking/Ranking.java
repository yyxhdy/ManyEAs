
package jmetal.util.ranking;

import jmetal.core.SolutionSet;

public interface Ranking {
	public SolutionSet getSubfront(int layer);
	public int getNumberOfSubfronts();
}
