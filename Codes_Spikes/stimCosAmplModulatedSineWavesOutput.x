#
#	stimCosAmplModulatedSineWavesOutput
#
#
foo nonperiodic	6
#
0	dur	"                Stimulus duration (s*10)"	20	1-600 1+
1	ampl	"               Max Amplitude (gain*1000)"	1000	0-1000 1+
2	SoundFreq	"            Frequency of pulses (kHz*10)"	1	1-260 1+
3	ModFreq	"  Amplitude Modulation Frequency (10*Hz)"	50	10-120 1+
4	Outputtype	"                 SoundCard(1), NIBord(2)"	2	1-2 1+
5	TriggerType	"     Manual(1), External(2),Immediate(3)"	3	1-3 1+
