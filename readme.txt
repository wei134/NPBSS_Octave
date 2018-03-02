====================================================================
		    
 NPBSS                              2018/03/01
 ver1.0.0                           ICI/LIFT, China
		    	 	       
====================================================================

System Requirements
===================

Runtime Environment Requirements
--------------------------------
Octave 4.2.1 or later


Software Requirements
---------------------
The NPBSS is supported on the following operating systems 
   
   Windows 7
   Windows 8 

Hardware Requirements/Recommendations
-------------------------------------
   Intel Core i3 2100 or later for Windows 
   Color display capable of 1024 X 768 pixel resolution
   

Memory Requirements/Recommendations
-------------------------------------
NPBSS requires approximately
5 G of available disk space on the drive or partition.

NPBSS requires a minimum of
2G of RAM for simulating PacBio CLR reads.


User's Guide
=================

Input Genome File Requirements
----------------------------
Input genome file is *.fasta or *.txt and input sequences require FASTA format, e.g.
 >Genome_name 
 AAACACTTGACATGTTCGTCGCGACTCTAAGAGATTAGAGTTTTCGGTTCGGCCGGACGAAACAC...


Execute Step
------------
Step 1: Open Octave.
Step 2: Change 'Current Folder' to the 'NPBSS' folder.
Step 3: Type in NPBSS('ecoilGenome.fa','-n 100') in the Command Window.
Step 4: Press 'ENTER' key to run.
Step 5: The final simulated sequences files can be found in the 'NPBSS' directory.

Test:
The scatter plot of read length and average base quality per read can be run by:
qv_average_each_read('simulated.fq')


Copyright Notice
===================
Software, documentation and related materials:
Copyright (c) 2018-2020 
Institute of Control & Information(ICI), Northwestern Polytechnical University,China
Key Laboratory of Information Fusion Technology(LIFT), Ministry of Education, China
All rights reserved.
