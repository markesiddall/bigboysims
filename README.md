# BIGBOYSIMS

## Description

This project was developed to examine the effects of various paramters on the 
ability to find the 'true' tree in the face of introgression/hybridization.
Specifically, it has been argued that finding the not-true tree is due to 
a history of introgression.  However, it is also possible that finding such a
not-true tree can be due to the use of a too-distant outgroup. 

This particular simulation uses the Drosophila melanogaster genome and its 
partitions into individual genes. In order to keep the size of the simulations
manageable, it uses the first 100 genes representing the first 240668 amino
acids. However, one could use all of them with two adjustments.
In def get_ungapped_sequence() you'll have to comment out this line
aligned_sequence = aligned_sequence[:240668]
And then under # parse the partitions:
change to file_path = 'melanogaster_partitions_all.txt'

### The simulated tree is:
```
    .----------... ? ...-------- OG
    |
    |    ,-------------- C
    |    |         |
    |    |         V
    `----|     ,-------- B
         `----|
               `-------- A
               
               <-depth->
          <----mya----->

mutation rate (mu) is varied
mya since MRCA is varied
depth of divergence of A and B is varied
the introgression --> is varied in amount from C to B, and occurs at half of depth
and the mutation mu*mya to the OG is varied according to mu and mya up to 200 Mya
```

### Execution to create simulations

If you run bigboy_simulation.py without an argument you will be prompted for how
many simulations you'd like to create and it will create them serially. 
Alternatively you can run bigboy_simulation.py with an argument and nest it in a 
loop to create a certian number each with a numerical index. This allows their
creation in parallel in the background...
for i in ${0..32}; do
    nohup bigboy_simulation.py $i &
    done; wait

### Results of simulation

Each simulation will result in a batch file for input to / execution by PAUP.
Given the suite of nested variables, there will be 1188000 simulated data sets.
Each data set has 14 taxa: A, B, C and 11 outgroups (OG). Iterating through 
partitions is accomlished by include/exclude commands in PAUP. Iterating over
OGs is accoplished by delete/restore commands in PAUP. 

### Execution in PAUP

You're on your own re PAUP execution. For example 
./paup4a169_ubuntu64 bigboy_sims_0.nex
nohup ./paup4a169_ubuntu64 bigboy_sims_0.nex &
The latter is recommended as it will take a while (~25 hr)
On AWS I'd recommend creating an AMI for PAUP, and then spin up as many
instances as you like (or are allowed) to run in parallel. 

For example
```
#!/bin/bash
for i in {0..32}
do
  aws ec2 run-instances \
    --image-id <your_id> \
    --instance-type c5.large \
    --key-name <your keypair> \
    --security-group-ids <your group id> \
    --subnet-id subnet-<your subnet>> \
    # here I had pushed the results of 32 simulations to an S3 bucket so they have to be pulled
    --user-data file://<(echo -e "#!/bin/bash\nN=$i\nsu - ubuntu -c 'aws s3 cp s3://bigboysims/results/bigboy_sims_$i.nex /home/ubuntu/bigboy_sims_$i.nex' \nsu - ubuntu -c '/home/ubuntu/paup4a169_ubuntu64 /home/ubuntu/bigboy_sims_$i.nex' \nsu - ubuntu -c 'aws s3 cp /home/ubuntu/simdump_$i.mrp s3://bigboysims/results/simdump_$i.mrp'") \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=PAUP_Instance_$i}]"
done
```

### Results of execution in PAUP

PAUP outputs the best 4-taxon tree in MRP format in a file named simdump_<index>.mrp, 
where <index> is whatever got passed to bigboy_simulation.py (e.g., numbers 0 through 31)
There will be 1188000 such mrps each appearing as such:
```
Begin data;
	Dimensions ntax=4 nchar=5;
	Format transpose;
	Taxlabels
		A
		B
		C
		OG2
		;
	Format transpose;
	Matrix

[Tree 1]
A       1000
B       0100
node_15 0011
C       0010
OG2     0001
	;
End;
```
If the line with 'node' is 0011 the tree is True
If the line with 'node' is 0110 the tree is Hybrid
If the line with 'node' is 0101 the tree is Wrong
If there is no line with 'node' the tree is Star

### Parsing MRP output from PAUP
Because the order of the MRPs is tied to the order of the simulations as such:
```
for mu in mu_values = [0.0001, 0.00025, 0.0005, 0.001, 0.0015, 0.002]
 for mya in mya_values = [2.5 * i for i in range(1, 11)]
   for depth in depth_values = [0.25, 0.5, 0.75]
    for proportion in introgression_proportions = [0, 0.01, 0.031, 0.0625, 0.125, 0.25]
     for partition in paritions 1 through 100
      for og in ogs = [2, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200]
```
the MRP output has to be faithful to that order

parse_mrps.py does this taking one simdump_<index>.mrp file as an argument. 
the output from simdump_<index>.mrp will be two files:

simdump_0.out contains all of the results for all partitions
simdump_0_summary.csv converts the results to have a count of number of partitions scoring

each of True, Hybrid, Wrong, and Star. 






