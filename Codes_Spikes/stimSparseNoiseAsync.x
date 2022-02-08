#
#	stimSparseNoiseAsync
#
#
foo nonperiodic	11
#
0	dur	"               Stimulus duration (s *10)"	10	1-60000 1+
1	x1	"       Left border, from center (deg*10)"	-1350	-1350-500 1+
2	x2	"      Right border, from center (deg*10)"	1350	-500-1350 1+
3	y1	"     Bottom border, from center (deg*10)"	-450	-1000-500 1+
4	y2	"        Top border, from center (deg*10)"	450	-500-1000 1+
5	nfr	"              Number of frames per image"	10	1-70 1+
6	sqsz	"                Size of squares (deg*10)"	30	1-200 1+
7	ncs	"  Number of contrasts (choose 2 for b&w)"	2	2-16 1+
8	c	"                            Contrast (%)"	100	0-100 1+
9	sprsns	"                Mean gray squares (%*10)"	950	0-1000 1+
10	seed	"         Seed of random number generator"	1	1-100 1+
