/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* Encapsulate a scoring matrix
*/

package edu.ucla.chem202.algorithm;

import edu.ucla.chem202.io.Logger;
import edu.ucla.chem202.DataFormatException;
	
import java.util.Arrays;
import java.util.List;

public class ScoringMatrix{
	
	private int[][] grid;
	private char[] headerList;	
	
	public ScoringMatrix(char[] header) {		
	}
	
	public void parseHeader(char[] header) {
		headerList = header;
		int headerSize = header.length;
		grid = new int[headerSize][headerSize];
	}
	
	private int getHeaderIndex(char toFind) {
		for (int count=0;count<headerList.length;count++) {
			if (toFind==headerList[count]) {
				return count;
			}
		}
		return -1;
	}
	
	public void parseLine(String[] line)
	throws DataFormatException{
		char toFind = line[0].charAt(0);
		int index = getHeaderIndex(toFind);
		if (index==-1) throw new DataFormatException ("The matrix residue headers do not match.  Be sure x and y axis headers are the same order.");
		if (line.length!=headerList.length+1) throw new DataFormatException("The number of data values in a row of the matrix is not correct.");
		for (int count=1;count<line.length;count++) {
			grid[index][count-1]=Integer.parseInt(line[count]);
		}
	}	
	
	public int getScore(char residue1, char residue2)
	throws DataFormatException{		
		int index1 = getHeaderIndex(residue1);
		int index2 = getHeaderIndex(residue2);
		if ((index1==-1) || (index2==-1)) throw new DataFormatException ("Residue in sequence not found in matrix.");
		//Logger.debugln("Trying to match a "+residue1+" to a "+residue2);
		return grid[getHeaderIndex(residue1)][getHeaderIndex(residue2)];
	}
}
