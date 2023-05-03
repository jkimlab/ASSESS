#!/bin/bash
if [ $# -lt 3 ] 
then
	echo "please check arguments"
	echo "bash ./ASSESS.MAKECMD.sh <input_dir> <threads> <output_dir> <running_assess>"
	echo "example: bash ./ASSESS.MAKECMD.sh ./ 10 test.out Y"
	exit
fi
INPUTDIR=$(realpath $1)
CPU=$2
OUTDIR=$3
RUN=$4
echo "## arguments"
echo "### input directory : $INPUTDIR"
echo "### num of threads : $CPU"
echo "### output directory : $OUTDIR"
mkdir $OUTDIR
echo "----------+"
date
echo "## [START] generating command lines using $INPUTDIR/params_runASSESSdocker.txt and $INPUTDIR/params_lib.txt"
docker run --rm -v $INPUTDIR:/assess_wd/run_test -t racaw:latest making_cmd.pl -lib /assess_wd/run_test/params_lib.txt -params /assess_wd/run_test/params_runASSESSdocker.txt -p $CPU > $OUTDIR/run_ASSESS_cmd.sh
echo "## [ END ] command lines are created [ $OUTDIR/run_ASSESS_cmd.sh ]"
echo "## [START] command line file post-processing"
sed -i 's/\r//g' $OUTDIR/run_ASSESS_cmd.sh
echo "## [ END ] all command lines are prepared"
if [ $# -eq 4 ]
then
	lcstr=$(echo $4 | tr '[:upper:]' '[:lower:]')
	if [ $lcstr = "y" ]
	then
		echo "### running option : $4"
		echo "## [START] running ASSESS"
		bash $OUTDIR/run_ASSESS_cmd.sh $OUTDIR > $OUTDIR/log.run_ASSESS.txt 2>&1
		echo "## [ END ] Final output : $OUTDIR/ASSESS_out"
	else
		echo "### running option : $4"
	fi
fi
date
echo "Finished"
