/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.gchen;

import edu.ucla.chem202.gchen.io.*;
import edu.ucla.chem202.gchen.algorithm.*;
import edu.ucla.chem202.gchen.algorithm.poa.*;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.List;
import java.util.Iterator;

public class ExerciseThree {
	
	private static void printUsage() {
		System.out.println("Usage: ./hbundle <clustal_infile> <perc_identity> <clustal_outfile>");
		System.exit(1);		
	}
	
	private void doConsensusGen(String args[])
	throws Exception{
		String clustalFileName = args[0];
		double percIdentity = Double.parseDouble(args[1]);
		String clustalOutFileName = args[2];
		Parser parser= new Parser();		
		List sequences = parser.parseClustal(new FileInputStream(clustalFileName));		
		if (sequences.size()<2) {
			Logger.println("Error: less than 2 sequences found, consensus generation aborting.");
			System.exit(1);
		}		
		Iterator it = sequences.iterator();
		while (it.hasNext()){
			Sequence seq = (Sequence)it.next();
		//	Logger.debugln("Sequence "+seq.getName()+" has contents "+seq.seqString());
		}
		POABundler poa = new POABundler(percIdentity);
		poa.build(sequences);
		FileOutputStream fos = new FileOutputStream(clustalOutFileName);
		poa.clustalOut(fos);
		fos.close();
		poa.outputBundleClassification();
	}

	public static void main(String args[]) {
		try {
			if (args.length<4) {
				printUsage();
			}
			int debugLevel = Integer.parseInt(args[3]);		
			Logger.setLevel(debugLevel);
			Logger.setOut(System.out);
			ExerciseThree exercise = new ExerciseThree();
			exercise.doConsensusGen(args);
		}
		catch (DataFormatException ex) {
			Logger.println("Formatting error: "+ex.getMessage());
		}
		catch (Exception ex) {
			ex.printStackTrace();
		}
	}	
}
