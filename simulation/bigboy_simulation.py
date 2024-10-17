
# %%

'''
The purpose of this simulation is to create batch input to PAUP so as to model introgression,
    and the relative success of finding the "true" tree under various parameter settings. 
    
The tree is:
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
and the mutation mu*mya to the OG is varied

'''
from sys import argv
import random
import copy

# Amino acids list (21 standard amino acids)
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Function to introduce mutations to an amino acid sequence
def mutate(sequence, mu, time):
    num_mutations = int(mu * time * len(sequence))  # Total mutations to apply
    sequence = list(sequence)  # Convert string to list for easy mutation
    
    for _ in range(num_mutations):
        pos = random.randint(0, len(sequence) - 1)
        sequence[pos] = random.choice(amino_acids)  # Mutate a random position with a random amino acid
    
    return ''.join(sequence)

# Function to apply introgression from C to B
def introgress(seq_from, seq_to, proportion):
    seq_from = list(seq_from)
    seq_to = list(seq_to)
    num_positions = int(proportion * len(seq_from))  # Number of positions to introgress
    
    positions = random.sample(range(len(seq_from)), num_positions)
    
    for pos in positions:
        seq_to[pos] = seq_from[pos]
    
    return ''.join(seq_to)

# Function to run the simulation for one set of parameters using the provided ungapped_sequence
def run_simulation_with_fixed_sequence(initial_seq, mu, mya, depth, introgression_proportion):
    B = initial_seq
    C = copy.deepcopy(B)  # Duplicate B to create C now (B C)
    
    # Mutate C for mu * mya
    C = mutate(C, mu, mya)
    
    # Mutate B for mu * (mya - depth), then split into A and B now ((A B) C)
    B = mutate(B, mu, mya*(1 - depth))
    A = copy.deepcopy(B)  # Create A from B
    A = mutate(A, mu, mya*depth)
    
    # Apply mutations to B for half of the depth before introgression
    half_depth = mya*depth/2
    B = mutate(B, mu, half_depth)
    
    # Apply introgression from C to B
    B = introgress(C, B, introgression_proportion)
    
    # Mutate B again for the remaining half of the depth
    B = mutate(B, mu, half_depth)
    
    return A, B, C

# Function to generate 11 additional sequences (OG2, OG5, etc.) mutated by different time values
def generate_OG_sequences(sequence, mu, times):
    og_sequences = {}
    for time in times:
        og_label = f"OG{time}"
        og_sequences[og_label] = mutate(sequence, mu, time)
    return og_sequences


#%% presently this does only one simulation

# get the partition names
def get_partition_names(file_path):
    partitions = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("charset"):
                # Extract the partition name (second string in the line)
                partition_name = line.split()[1]
                partitions.append(partition_name)
    
    return partitions

# get the genome
def get_ungapped_sequence(file_path):
    with open(file_path, 'r') as f:
        # Read the first line (the aligned sequence)
        aligned_sequence = f.readline().strip()
        
        # Since this is already ungapped (except '-') just return it
        aligned_sequence = aligned_sequence[:240668]
        return aligned_sequence


# Function to output each simulation result to the NEXUS file
def output_simulation_to_nexus(simulation_results, paupfile, idx, mu, mya, depth, introgression_proportion,sim):
    
    # create treefile name
    treefile = 'simdump_'+str(sim)+'.mrp'
    # Write simulation info
    print(f"[Simulation {idx+1}: mu= {mu}, mya= {mya}, depth= {depth}, introgression_proportion= {introgression_proportion}]", file=paupfile)
    print("Begin data;", file=paupfile)
    print("Dimensions ntax=14 nchar=", nchar ,";", file=paupfile)
    print("Format datatype=protein gap=-;", file=paupfile)
    print("Matrix", file=paupfile)
    
    # Write the sequences A, B, C
    print(f"A {simulation_results['A']}", file=paupfile)
    print(f"B {simulation_results['B']}", file=paupfile)
    print(f"C {simulation_results['C']}", file=paupfile)
    
    # Write OG sequences (mutated from the original sequence at different times)
    for og_label, og_sequence in simulation_results['OG_sequences'].items():
        print(f"{og_label} {og_sequence}", file=paupfile)
    
    print("; end;", file=paupfile)
    
    for partition in partitions:
        print("Begin paup;", file=paupfile)
        print("execute melanogaster_partitions.txt;", file=paupfile)
        print("exclude all; include", partition, ";", file=paupfile)
        for og in range (4,15):
            print("delete all /cleartrees=yes; restore 1-3;  restore", og,";", file=paupfile)
            print("alltrees; matrixrep file=", treefile, "append=yes;", file=paupfile)
        print("end;", file=paupfile)
    

# Function to loop through simulations and write results to the NEXUS file
def run_and_save_simulations_nexus(mu_values, mya_values, depth_values, introgression_proportions, output_file,sim):
    with open(output_file, "w") as paupfile:
        # Write NEXUS preamble
        print("#NEXUS", file=paupfile)
        print("begin paup;", file=paupfile)
        #print("cd /Users/msiddall/Documents/RobTriplets/simulations/;", file=paupfile)
        print("set monitor = no;", file=paupfile)
        print("set warnReset = no;", file=paupfile)
        print("set warnTree  = no;", file=paupfile)
        print("end;", file=paupfile)
        
        idx = 0
        
        for mu in mu_values:
            for mya in mya_values:
                for depth in depth_values:
                    for introgression_proportion in introgression_proportions:
                        # Increment simulation counter
                        idx += 1
                        
                        # Run the simulation
                        A, B, C = run_simulation_with_fixed_sequence(ungapped_sequence, mu, mya, depth, introgression_proportion)
                        
                        # Generate additional OG sequences (mutated by different times)
                        og_times = [2, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200]
                        og_sequences = generate_OG_sequences(ungapped_sequence, mu, og_times)
                        
                        # Store the simulation results
                        simulation_results = {
                            "A": A,
                            "B": B,
                            "C": C,
                            "OG_sequences": og_sequences
                        }
                        
                        # Write the results to the NEXUS file
                        output_simulation_to_nexus(simulation_results, paupfile, idx, mu, mya, depth, introgression_proportion,sim)
                paupfile.flush()

# The suite of parameters (save for outgroups)
mu_values = [0.0001, 0.00025, 0.0005, 0.001, 0.0015, 0.002]
mya_values = [2.5 * i for i in range(1, 11)]
depth_values = [0.25, 0.5, 0.75]
introgression_proportions = [0, 0.01, 0.031, 0.0625, 0.125, 0.25]


# parse the genome
file_path = 'melanogaster_genome.txt'
ungapped_sequence = get_ungapped_sequence(file_path)
nchar = len(ungapped_sequence)
# parse the partitions:
file_path = 'melanogaster_partitions.txt'
partitions = get_partition_names(file_path)

# Run and save simulations to NEXUS file

# in the case of interactive  and serial
if len(argv) == 0:
    num_sims=input("Number of simulations: ")
    num_sims=int(num_sims)
    for sim in range (0, num_sims):
        output_file = 'bigboy_sims_'+str(sim)+'.nex'
        run_and_save_simulations_nexus(mu_values, mya_values, depth_values, introgression_proportions, output_file, sim)

# in the case of an external bash or Nextflow loop passing in parallel: 
#      for i in {0..99}
#         do
#           python3 bigboy_simulation.py $i &
#         done
else:
    sim=int(argv[1])
    output_file = 'bigboy_sims_'+str(sim)+'.nex'
    run_and_save_simulations_nexus(mu_values, mya_values, depth_values, introgression_proportions, output_file, sim)

