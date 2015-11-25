#!/usr/bin/perl

our %table;


while (<STDIN>){

    if(/[0-9\.-]+\s+\S+\s+(\S+)\s+/){
        #print "$1\n";
        $table{$1}=0;
    }
    
}

# LOOP THROUGH IT
while (($key, $value) = each(%table)){
     #print $key."  ".$value."  \n";
    print "      <struct name=\"FieldValueShape$key\">\n";
    print "         <param name=\"Type\">FieldValueShape</param>\n";
    print "         <param name=\"ValueField\">PlateID_VoxelField</param>\n";
    $keyL=$key-0.5;
    $keyU=$key+0.5;
    print "         <param name=\"LowerLimit\">$keyL</param>\n";
    print "         <param name=\"UpperLimit\">$keyU</param>\n";
    print "      </struct>\n";

    print "      <struct name=\"viscosity$key\">\n";
    print "         <param name=\"Type\">MaterialViscosity</param>\n";
    print "         <param name=\"eta\">$key</param>\n";
    print "      </struct>\n";

    print "      <struct name=\"material$key\">\n";
    print "         <param name=\"Type\">RheologyMaterial</param>\n";
    print "         <param name=\"Shape\">FieldValueShape$key</param>\n";
    print "         <param name=\"density\">1.0</param>\n";
    print "         <param name=\"Rheology\">viscosity$key</param>\n";
    print "      </struct>\n";

}
#while (($key, $value) = each(%table)){
#    print $key."  ".$value."  \n";
#}

