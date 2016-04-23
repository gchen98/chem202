/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.gchen;

import edu.ucla.chem202.gchen.io.*;
import edu.ucla.chem202.gchen.algorithm.*;
import edu.ucla.chem202.gchen.algorithm.dp.*;
import java.io.FileInputStream;
import java.util.List;
import java.util.Iterator;

public class ExerciseTwo {
	
	private static void printUsage() {
		System.out.println("Usage: java ExerciseTwo <matrix filename> <ealign matrix> "+
			"<fasta sequence filenam> <gap open> <gap extension> <debug_level>");
		System.exit(1);		
	}
	
	private void doGlobal(String args[])
	throws Exception{
		String matrixFileName = args[0];
		String fastaFileName = args[2];
		int gapOpening = Integer.parseInt(args[3]);
		int gapExtend = Integer.parseInt(args[4]);		
		Parser parser= new Parser();		
		ScoringMatrix matrix = parser.parseMatrix(new FileInputStream(matrixFileName));		
		List sequences = parser.parseFasta(new FileInputStream(fastaFileName));
		if (sequences.size()>2) Logger.debugln("Warning: more than 2 sequences found, only first two will be used.");
		if (sequences.size()<2) {
			Logger.println("Error: less than 2 sequences found, alignment aborting.");
			System.exit(1);
		}		
		Aligner aligner= new TwoSequenceGlobalAligner();
		aligner.init(sequences,	gapOpening,gapExtend, matrix);
		aligner.align();			
		aligner.dump();
		Alignment alignment = aligner.getAlignment();
		Logger.println("CLUSTAL format: ");
		alignment.toClustal();
		Logger.println("Global alignment has score: "+aligner.getScore());		
	}
	
	private void doHairpin(String args[])
	throws Exception{
		String invertMatrixFileName = args[1];
		String fastaFileName = args[2];
		int gapOpening = Integer.parseInt(args[3]);
		int gapExtend = Integer.parseInt(args[4]);
		Parser parser= new Parser();
		List sequences = parser.parseFasta(new FileInputStream(fastaFileName));
		if (sequences.size()>1) Logger.debugln("Warning: more than 1 sequences found, only first one will be used.");
		if (sequences.size()<1) {
			Logger.println("Error: less than 1 sequences found, alignment aborting.");
			System.exit(1);
		}
		Sequence einvertSeq = (Sequence)sequences.get(0);
		Sequence revCompl = new Sequence(einvertSeq.getName(),einvertSeq.getDesc(),
			einvertSeq.reverseComplement());
		sequences.set(1,revCompl);
		ScoringMatrix matrix = parser.parseMatrix(new FileInputStream(invertMatrixFileName));
		//matrix = parser.parseMatrix(new FileInputStream(matrixFileName));
		Aligner aligner= new TwoSequenceLoopAligner();
		aligner.init(sequences, gapOpening,gapExtend, matrix);
		aligner.align();
		aligner.dump();
		Alignment alignment = aligner.getAlignment();
		Logger.println("CLUSTAL format: ");
		alignment.toClustal();
		Logger.println("Hairpin alignment has score: "+aligner.getScore());		
	}
	
	public static void main(String args[]) {
		try {
			if (args.length<6) {
				printUsage();
			}
			int debugLevel = Integer.parseInt(args[5]);		
			Logger.setLevel(debugLevel);
			Logger.setOut(System.out);
			ExerciseTwo exercise = new ExerciseTwo();
			//exercise.doGlobal(args);
			exercise.doHairpin(args);
		}
		catch (DataFormatException ex) {
			Logger.println("Formatting error: "+ex.getMessage());
		}
		catch (Exception ex) {
			ex.printStackTrace();
		}
	}	
}
