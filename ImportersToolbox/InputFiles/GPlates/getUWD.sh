str=`pwd`

arr1=( `echo "$str" | tr -s '/' ' '` )
UWDIR=""

#for name in ${arr1[@]}
#do
#echo $name
# other stuff on $name
#done

for (( i = 0 ; i < ${#arr1[@]} ; i++ ))
do
if test  ${arr1[$i]} ==  "ImportersToolbox"
then
    #echo "Found Underworld at level $i"
    #construct path to UW
    for (( j = 0 ; j < $i ; j++ ))
    do
	UWDIR=$UWDIR/${arr1[$j]}
	#echo "$UWDIR"
    done
#else
    #echo "Underworld Not Found"
fi
done

echo "$UWDIR"
#echo "from getUWD -- $0"