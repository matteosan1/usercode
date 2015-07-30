#! /bin/bash

# send.sh queue directory_name n_events start end

echo 'Sending jobs...'

for i in `seq $4 $5`;
do
 bsub -q $1 job1.sh $2 $3 $i
done    
