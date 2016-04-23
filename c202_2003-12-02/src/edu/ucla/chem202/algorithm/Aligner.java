/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.algorithm;

import edu.ucla.chem202.*;
import java.util.List;

/**
* Align two sequences using dynamic programming.
*/
	
public interface Aligner{
	
	// Let's assume sequence 1 to rest on the y axis and sequence 2 on the x axis
	
	public void init(List sequences, int gapOpen, int gapExtend, ScoringMatrix matrix);
	
	public int getScore();	
	
	public Alignment getAlignment();
	
	public void align()
	throws DataFormatException;
	
	public void dump();	
}
