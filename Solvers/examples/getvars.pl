#!/usr/bin/perl

$str=$ARGV[0];

while(<STDIN>){

    if(/$str/){
	s/\|//g;
	print $_;
    }

}
