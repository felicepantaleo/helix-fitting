#!/usr/local/bin/perl

system ("rm data/*");
$MANY=2000;
#$MANY=1;
# @pTs=("0_8", "1_5", "2_5", "5", "10", "20", "50" );
@pTs=("1", "1_5", "2_5", "5", "10", "20", "50" );
#@pTs=("1");
@types=("barrel","forward" );
@fits=("R0", "G0", "K0", "R1", "G1", "K1", "R2" , "G2", "G3");
# @pTs=("0_8");
@fitlong=("riemann0", "global0", "kalman0", "riemann1", "global1", "kalman1" , "riemann2", "global2", "global3" );
$seed=7269336;
$cm="traxx -n $MANY ";
	# $type="barrel";
foreach $type (@types) {
	foreach $pT (@pTs) {
		foreach $ft (@fits) {
			$command=$cm."-".$ft." -f etc/$type.conf -i ".
				"etc/$type$pT.init.conf -s$seed --delphi=data/$type$pT -t";
			print "issuing: $command\n";
			system($command);
		};
	};
};
