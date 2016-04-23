package edu.ucla.chem202.io;

import java.io.PrintStream;

/*
 * Author: gchen
 * Created: Wednesday, October 15, 2003 7:17:52 AM
 * Modified: Wednesday, October 15, 2003 7:17:52 AM
 */


public class Logger
{
	public static final int QUIET=0;
	public static final int DEBUG=1;
	
	private static PrintStream out;
	private static int level;
	
	public static void setLevel(int level) {
		Logger.level = level;
	}
	
	public static int getLevel() {
		return level;
	}
	
	
	public static void setOut(PrintStream out) {
		Logger.out = out;
	}	

	private static void output(String str, int level, boolean newLine){
		if (level<=Logger.level) {
			if (newLine) out.println(str);
			else out.print(str);
		}
	}
	
	public static void debug(String str) {
		output(str,DEBUG,false);
	}
	
	public static void debugln(String str){
		output(str,DEBUG,true);
	}	
		
	public static void print(String str) {
		output(str,QUIET,false);
	}
	
	public static void println(String str) {
		output(str,QUIET,true);
	}	
	
}
