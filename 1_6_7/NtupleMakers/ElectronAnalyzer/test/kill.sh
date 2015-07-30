#! /bin/bash

for i in `seq $1 $2`;
do
 echo 'Killing job n. $i'
 bkill $i
done    
