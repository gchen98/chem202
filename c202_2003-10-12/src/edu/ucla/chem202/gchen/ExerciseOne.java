/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.gchen ;

import edu.ucla.chem202.gchen.io.*;
import edu.ucla.chem202.gchen.algorithm.*;
import edu.ucla.chem202.gchen.algorithm.dp.*;

import java.io.FileInputStream;
import java.util.List;
import java.util.Iterator;

public class ExerciseOne {
	
	private static void printUsage() {
		System.out.println("Usage: java ExerciseOne <matrix filename> "+
			"<fasta sequence filenam> <gap open> <gap extension> <debug_level>");
		System.exit(1);		
	}
	
	public static void main(String args[]) {
		try {	
			if (args.length<5) {
				printUsage();
			}
			String matrixFileName = args[0];
			String fastaFileName = args[1];
			int gapOpening = Integer.parseInt(args[2]);
			int gapExtend = Integer.parseInt(args[3]);
			int debugLevel = Integer.parseInt(args[4]);
			
			Logger.setLevel(debugLevel);
			Logger.setOut(System.out);
			
			Parser parser= new Parser();
			
			ScoringMatrix matrix = parser.parseMatrix(new FileInputStream(matrixFileName));
			
			List sequences = parser.parseFasta(new FileInputStream(fastaFileName));
			if (sequences.size()>2) Logger.println("Warning: more than 2 sequences found, only first two will be used.");
			if (sequences.size()<2) {
				Logger.println("Error: less than 2 sequences found, alignment aborting.");
				System.exit(1);
			}
			
			Aligner aligner= new TwoSequenceLocalAligner();
			//Aligner aligner= new ObjectCellLocalAligner();
			aligner.init(sequences, gapOpening,gapExtend, matrix);
			aligner.align();
			aligner.dump();
			Alignment alignment = aligner.getAlignment();			
			Logger.println("CLUSTAL format: ");
			alignment.toClustal();
			Logger.println("Alignment has score: "+aligner.getScore());
		}
		catch (DataFormatException ex) {
			Logger.println("Formatting error: "+ex.getMessage());
		}
		catch (Exception ex) {
			ex.printStackTrace();
		}
	}	
}
