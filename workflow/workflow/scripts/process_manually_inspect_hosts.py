#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---

##################################################
# IMPORTS
##################################################
from datetime import datetime


##################################################
# FUNCTIONS
##################################################
def process_mappings(mapping_path):
    import pandas as pd
    import numpy as np

    special_category = ["unknown", "environmental"]

    mapping = pd.read_csv(mapping_path)
    # Complete Eval column
    mapping.loc[:, "Eval"] = mapping["Eval"].astype("boolean")
    mapping.loc[:,"Eval"] = np.where(mapping["Eval"] == False,
                                     False, True)
    # Complete Final_host column

    mapping.loc[:, "Final_host"] = np.where(mapping["Eval"] == True,
                                    np.where(mapping['Assigned_method']=="environmental",
                                            "environmental",
                                            np.where(mapping['Final_host'].isin(special_category),
                                                mapping["Final_host"],
                                                mapping["Assigned_host"])),
                                    mapping["Final_host"])
    mapping.to_csv(mapping_path, index=False)
    
    return mapping
##################################################
# ARGS
##################################################
def get_parser():
    import argparse
    # create parser object
    parser = argparse.ArgumentParser(
        description="Ensure that manually inspection of mappings has been performed",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # specify our desired options 
    # by default ArgumentParser will add an help option 
    parser.add_argument("--input", help="csv file for host mapping from BIOSAMPLE_Host")
    parser.add_argument("--outfile")

    return parser

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    ARGS = get_parser().parse_args()

    # Define the task
    task_name = "Manually inspection of mapping file"
    task = f"""This is the task: Manually inspection of the mapping file:
        - {ARGS.input}
         1) Select the new entries of the mapping (filter by version)
         2) Select the entries different from the Assigned_methods "ncbi_taxonomy" and "ete3" (we assume they are correct)
         3) Evaluate if the Assigned_host is correct:
            - If yes, do nothing
            - If no, write False in the 'Eval' column and write the correct host name in 'Final_host' column
            - If 'Assigned_method' == 'not_assigned':
               - and you think the sample is environmental, write environmetal in the 'Final_host' column
               - If you don't find an appropiate host, write 'unknown' in 'Final_host' column
         4) Please, add new rules to 'fix_host' (process_create_mappings.py) to improve the assignment (if convinient)
         5) Save the document and continue the pipeline :D
                """

    # Ask the user if the task has been completed
    print(task)
    user_input = input(f"Have you completed the {task_name} ? (yes/no): ").strip().lower()

    # Check the user's response
    if user_input == "yes":
        print(f"Great job! You've completed '{task_name}'.")
        
        # Process mappings
        mapping = process_mappings(ARGS.input) 
        
        # save mappings
        # Record datetime dd/mm/YY H:M:S
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        # Create outputfile
        with open(ARGS.outfile, "w+") as out:
            out.write(f"[{dt_string}] The following task has been completed :\n")
            out.write(task)
        
    elif user_input == "no":
        print(f"Don't worry, you can still complete '{task_name}'. Come back again after completion")
    else:
        print("Invalid response. Please enter 'yes' or 'no'.")