/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202;

import edu.ucla.chem202.io.*;
import edu.ucla.chem202.algorithm.*;
import edu.ucla.chem202.algorithm.hmm.*;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.List;
import java.util.Iterator;

public class ExerciseFive {
	
	private static void printUsage() {
		System.out.println("Usage: ./profileprob <clustal_infile> <fasta_infile> <residue_definition_file> <debug_level(0=terse,1=verbose)>");
		System.exit(1);		
	}
	
	private void execute(String clustalFileName, String fastaFileName, String residueDefFileName)
	throws Exception{
		Parser parser= new Parser();
		List clustalSequences = parser.parseClustal(new FileInputStream(clustalFileName));
		List fastaSequences = parser.parseFasta(new FileInputStream(fastaFileName));
		HMM hmm = new HMMAlignment(clustalSequences, new FileInputStream(residueDefFileName));
		HMMAligner aligner = new HMMGlobalAligner(); // replace null with the name of Teresa's new class
		aligner.init(hmm);
		Iterator it = fastaSequences.iterator();
		while (it.hasNext()) {
			Sequence seq = (Sequence)it.next();
			aligner.align(seq);
			Logger.println(seq.getName()+" LOG_PROB="+aligner.getTotalProb());
			if (Logger.getLevel()>=Logger.DEBUG) {
				aligner.printDataStructure();
			}
		}
	}

	public static void main(String args[]) {
		try {
			if (args.length<4) {
				printUsage();
			}
			int debugLevel = Integer.parseInt(args[args.length-1]);			
			Logger.setLevel(debugLevel);
			Logger.setOut(System.out);
			ExerciseFive exercise = new ExerciseFive();
			exercise.execute(args[0], args[1], args[2]);
		}		
		catch (Exception ex) {			
			Logger.println("Application error of type "+ex.getClass().getName()+": "+ex.getMessage());
			if (Logger.getLevel()>=Logger.DEBUG) {
				ex.printStackTrace(System.out);
			}
		}
	}	
}
