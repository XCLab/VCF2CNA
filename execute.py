#!/usr/bin/python

# To run this module, you must make sure python 2.7 is running

import argparse
import sys, os
import subprocess

def parseArguments():
	# Create argument parser
	parser = argparse.ArgumentParser()

	# mandatory argument
	parser.add_argument("FileName", help="File to Analyze.",  type=str)
	parser.add_argument("FilePath", help="Path to FileName.", type=str)
	parser.add_argument("OutputDirectory" , help="Path to Store files.", type=str)


	# optional argument (The "-" prefix of the argument name specifies optional status)
	parser.add_argument("-VcfOrder", help="Order of paired VCF File (TN or NT).", type=str, default="TN")
	parser.add_argument("-DiploidChr", help="Chromosome used to normalize results", type=str, default="-1")
	parser.add_argument("-StartChr", help="Start of chr segment on diploid chr to normalize results", type=str, default="-1")
	parser.add_argument("-EndChr", help="End of chr segment on diploid chr to normalize results", type=str, default="-1")
	parser.add_argument("-Median", help="Median Normal Coverage", type=str, default="-1")
	parser.add_argument("-MinSf", help="Minimum Scale Factor", type=str, default="0.50")
	parser.add_argument("-MaxSf", help="Maximum Scale Factor", type=str, default="1.50")
	parser.add_argument("-XminSf", help="Minimum X Scale Factor", type=str, default="0.25")
	parser.add_argument("-XmaxSf", help="Maximum X Scale Factor", type=str, default="1.50")
	parser.add_argument("-R_Loc",    help="Path to R Program", type=str, default="")
	parser.add_argument("-R_ScriptLoc", help="Path to RScript Program", type=str, default="")


	# Parse arguments
	args = parser.parse_args()
	return args

def main(LineArgs):
	filename = LineArgs.FileName
	pathname = LineArgs.FilePath
	workdir  = LineArgs.OutputDirectory
	R_Loc    = LineArgs.R_Loc
	R_ScriptLoc = LineArgs.R_ScriptLoc
	VcfOrder = LineArgs.VcfOrder
	DiploidChr = LineArgs.DiploidChr
	Median     = LineArgs.Median 
	MinSf      = LineArgs.MinSf
	MaxSf      = LineArgs.MaxSf
	XminSf     = LineArgs.XminSf
	XmaxSf     = LineArgs.XmaxSf
	StartChr   = LineArgs.StartChr
	EndChr     = LineArgs.EndChr

	args=['./commandline.sh', filename, pathname, VcfOrder, DiploidChr, Median, MinSf, MaxSf, XminSf, XmaxSf, workdir, R_Loc, R_ScriptLoc, StartChr, EndChr]

	p = subprocess.Popen(args)
	p.wait()
	exit()

if __name__ == '__main__':
	# Parse args
	LineArgs = parseArguments()
	main(LineArgs)

