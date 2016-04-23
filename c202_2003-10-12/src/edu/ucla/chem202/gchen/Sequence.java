/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* Encapsulate a sequence and its attributes
*/

package edu.ucla.chem202.gchen;
	
import java.io.InputStream;
import java.util.HashMap;

public class Sequence{
	
	private String name;
	private String desc;
	private String seq;
	private double identityToCons;
	public static HashMap bases = new HashMap();
	
	public Sequence(String name, String desc, String seq) {
		this.name = name;
		this.desc = desc;
		seq.replace('_','-');		
		this.seq = seq.toUpperCase();
		bases.put("A","T");
		bases.put("T","A");
		bases.put("C","G");
		bases.put("G","C");
		bases.put("N","N");
	}
	
	
	public String getName(){
		return name;
	}
	
	public void setIdentityToConsensus(double perc) {
		this.identityToCons = perc;
	}
	
	public double getIdentityToConsensus() {
		return this.identityToCons;
	}
	
	public String getDesc() {
		return desc;
	}
	
	public void setDesc(String desc) {
		this.desc = desc;
	}
	
	public String seqString() {
		return seq;
	}
	
	public String complement() {
		StringBuffer buf = new StringBuffer();
		for (int i=0;i<seq.length();i++) {
			buf.append(bases.get(seq.substring(i,i+1)).toString());
		}
		return buf.toString();
	}
	
	public String reverseComplement() {
		StringBuffer buf = new StringBuffer();
		for (int i=seq.length()-1;i>=0;i--) {
			buf.append(bases.get(seq.substring(i,i+1)).toString());
		}
		return buf.toString();
	}
		
	public static void main(String args[]) {
		Sequence sequence = new Sequence("gary","desc","AAACCCTTTGGTTTTTT");
		System.out.println(sequence.complement());
		System.out.println(sequence.reverseComplement());
	}

}
