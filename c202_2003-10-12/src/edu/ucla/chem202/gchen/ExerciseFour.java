/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.gchen;

import edu.ucla.chem202.gchen.io.*;
import edu.ucla.chem202.gchen.algorithm.*;
import edu.ucla.chem202.gchen.algorithm.poa.*;
import edu.ucla.chem202.gchen.algorithm.dp.*;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.List;
import java.util.Iterator;

public class ExerciseFour {
	
	private static void printUsage() {
		System.out.println("Usage: ./poalign <fasta_infile> <matrix_file> <gap_penalty> <clustal_outfile> <debug_level(0=terse,1=verbose)>");
		System.exit(1);		
	}
	
	private void doAlign(String clustalFileName, String matrixFileName,
		int gapPenalty, String outputFileName)
	throws Exception{
		Parser parser= new Parser();
		ScoringMatrix matrix = parser.parseMatrix(new FileInputStream(matrixFileName));
		List sequences = parser.parseFasta(new FileInputStream(clustalFileName));
		POALocalAligner aligner = new POALocalAligner();
		aligner.init(sequences, gapPenalty, gapPenalty, matrix);
		aligner.align();
		aligner.clustalOut(new FileOutputStream(outputFileName));
	}

	public static void main(String args[]) {
		try {
			if (args.length<5) {
				printUsage();
			}
			int debugLevel = Integer.parseInt(args[4]);
			Logger.setLevel(debugLevel);
			Logger.setOut(System.out);
			ExerciseFour exercise = new ExerciseFour();
			exercise.doAlign(args[0], args[1],Integer.parseInt(args[2]), args[3]);
		}
		catch (DataFormatException ex) {
			Logger.println("Formatting error: "+ex.getMessage());
		}
		catch (Exception ex) {
			ex.printStackTrace();
		}
	}	
}
