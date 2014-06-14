#!/usr/bin/perl
$sump=0;
$sumv=0;
$timep=0;
$nlits=0;
  
while(<STDIN>){
    if(/Pressure Solve:         = (\d.*\d) secs \/ (\d+) its/){
        $sump = $sump + $2;
	$timep = $timep + $1;
	$nlits = $nlits+1;
    }
    if(/Final V Solve:\s+= (\d.*\d) secs \/ (\d+) its/){
	#$sumv = $sumv + $2;
    }
    if(/Linear solve converged due to CONVERGED_RTOL iterations (\d+)/){
	$sumv = $sumv + $1;
    }
    if(/min\/max\(u\)    = (\-\d.*\d) \[\d+\] \/ (\d.*\d) \[\d+\]/){
	$minu=$1; $maxu=$2;
    }
    if(/min\/max\(p\)    = (\d.*\d) \[\d+\] \/ (\d.*\d) \[\d+\]/){
	$minp=$1; $maxp=$2;
    }
    if(/min\/max\(p\)    = (\-\d.*\d) \[\d+\] \/ (\d.*\d) \[\d+\]/){
	$minp=$1; $maxp=$2;
    }
    if(/PressureField - dof 0 global error: (\d.*\d)/){
	$gperr=$1;
    }
    if(/= (\d.*\d)\s+:constraint/){
	$pdoferr=$1;
    }
}
#print "Total Pressure its: $sump  Total P time = $timep\n";  
#print "Total Velocity its: $sumv\n";  
#print "min U: $minu : max U: $maxu\n";
#print "min P: $minp : max P: $maxp\n";
#|G^T u + C p - h|        = 2.62589986e-03  :constraint
$avp=$sump/$nlits;
print "| $sump | $sumv | $timep | $pdoferr | $gperr | $nlits | $avp |$minp | $maxp | $minu | $maxu |"
