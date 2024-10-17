# %%

import sys
import re
import pandas as pd

#%%
if __name__ == '__main__':
    if len(sys.argv) < 3:
        filename = input ("Enter input file: ")
    else:
        filename = sys.argv[1]
    outfilename = filename[0:-4] + ".out"  
    # %%
    # Parameter values
    mu_values = [0.0001, 0.00025, 0.0005, 0.001, 0.0015, 0.002]
    mya_values = [2.5 * i for i in range(1, 11)]
    depth_values = [0.25, 0.5, 0.75]
    introgression_proportions = [0, 0.01, 0.031, 0.0625, 0.125, 0.25]
    partitions = range(100)
    og_times = [2, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200]
    
    node_pattern = re.compile(r"node[^a-zA-Z]")

    # Open the input file for reading and output file for writing
    with open(filename, 'r') as file, open(outfilename, 'w') as outfile:
        # Write the header to the output file
        outfile.write("mu,mya,depth,introgression_proportion,partition,og_time,tree_type,True,Hybrid,Wrong,Star\n")
        
        tree_index = 0
        
        # Iterate over each parameter in the specified order
        for mu in mu_values:
            for mya in mya_values:
                for depth in depth_values:
                    for introgression_proportion in introgression_proportions:
                        for partition in partitions:
                            for og_time in og_times:
                                # Read the next tree block from the file
                                tree_block = []
                                while True:
                                    line = file.readline()
                                    if not line or line.startswith('End;'):
                                        break
                                    tree_block.append(line)
                                
                                # Determine if the tree is T, H, W, or S
                                tree_type = "S"  # Default to 'S' if no node is found
                                for line in tree_block:
                                    if node_pattern.search(line):
                                        if "0011" in line:
                                            tree_type = "T"
                                        elif "0110" in line:
                                            tree_type = "H"
                                        elif "0101" in line:
                                            tree_type = "W"
                                        break  # We can stop processing once we've determined the type
                                
                                # Assign the values for the new columns
                                is_true = 1 if tree_type == "T" else 0
                                is_hybrid = 1 if tree_type == "H" else 0
                                is_wrong = 1 if tree_type == "W" else 0
                                is_star = 1 if tree_type == "S" else 0

                                # Write result to output file
                                outfile.write(f"{mu},{mya},{depth},{introgression_proportion},{partition},{og_time},{tree_type},{is_true},{is_hybrid},{is_wrong},{is_star}\n")
                                
                                tree_index += 1  # Move to the next tree

    # Output completion message
    print(f"Results written to {outfilename}")
    
    # File paths
    input_filename = outfilename
    output_summary_filename = outfilename[:-4]+"_summary.csv"
    
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_filename)
    
    # Group by the parameters of interest and sum the "True", "Hybrid", "Wrong", and "Star" columns
    summary_df = df.groupby(['mu', 'mya', 'depth', 'introgression_proportion', 'og_time']).agg(
        True_count=('True', 'sum'),
        Hybrid_count=('Hybrid', 'sum'),
        Wrong_count=('Wrong', 'sum'),
        Star_count=('Star', 'sum')
    ).reset_index()
    
    # Write the summary DataFrame to a new CSV file
    summary_df.to_csv(output_summary_filename, index=False)

    # Output completion message
    print(f"Summary results written to {output_summary_filename}")

