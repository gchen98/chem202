#!/usr/bin/perl -w
use strict;

if (@ARGV<1) {
 print "Enter debug level\n";
 exit(1); 
}

#global constants
my $inputdir = "./input_sequences";
my $debuglevel = shift;

sub testCase{
 my $comment = shift;
 my $clustalFile = shift;
 my $fastaFile = shift;
 my $residueFile = shift;
 print "\nRunning test case for $comment\n****************************\n";
 my $command = "./posteriorprob $inputdir/$clustalFile $inputdir/$fastaFile $inputdir/$residueFile $debuglevel";
 system($command);
 print "****************************\nEnd test case for $comment\n\n";
}

#testCase("normal case med seq","mediumAlign.txt","mediumSeq.txt","amino_acids.txt");
#testCase("normal case long dna","alignment.txt","sequences.txt","nucleotides.txt");
#testCase("bad alignment file","asyAlign.txt","easySeq.txt","amino_acids.txt");
#testCase("all same chars","mediumAlign.txt","sameCharSeq.txt","amino_acids.txt");
#testCase("bad residue in seq","mediumAlign.txt","mediumSeq.txt","nucleotides.txt");
#testCase("long-assed sequence","globinAlign.txt","globinSeq1.txt","amino_acids.txt");
#testCase("should insert in beginning of alignment","teresaAlign.txt","teresaSeq.txt","amino_acids.txt");
#testCase("adeno alignment","Globin.txt","mikeSeq.txt","amino_acids.txt");
#testCase("simple case short seq","easyAlign.txt","easySeq.txt","amino_acids.txt");
#testCase("medium length alignment","mediumAlign.txt","mediumSeq.txt","amino_acids.txt");
#testCase("mike's alignment","mikeAlign.txt","mikeSeq.txt","amino_acids.txt");
testCase("stuart's","stuartAlign.txt","stuartSeq.txt","amino_acids.txt");
