/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.algorithm.hmm;

import edu.ucla.chem202.*;
import edu.ucla.chem202.algorithm.Sequence;
import java.util.List;

/**
* Align two sequences using dynamic programming.
*/
	
public interface HMMAligner{
	
	public void init(HMM hiddenMarkovModel);	
	
	public void align(Sequence seq)
	throws DataFormatException;
	
	public double getTotalProb();
	
	public void printDataStructure();
	
	public void printPosteriorProbs();
}