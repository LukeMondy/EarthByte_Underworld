#!/usr/bin/perl
$sum=0;
print "\n";                                                                                                                  
while(<STDIN>){
    if(/Non linear solver - Residual (\d.*\d); Tolerance (\d.*\d) - Converged - (\d.*\d) \(secs\)/){
	print "Solution time:      $3\n";
	print "Picard residual:    $1\n";
    }
    if(/In func SystemLinearEquations_NonLinearExecute: Converged after (\d+) iterations/){
	print "Picard Iterations:  $1\n";
    }
    if(/Pressure Solve:         = (\d.*\d) secs \/ (\d+) its/){
        $sum = $sum + $2;
    }                                                                                                                                   
    if(/Pressure Solve P its =\s+(\d+)/){
	$sum = $sum + $1;
    }
}
print "Total Pressure its: $sum\n";  
