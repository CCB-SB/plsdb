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
    parser.add_argument("--input-locations", help="csv file for host mapping from BIOSAMPLE_Host")
    parser.add_argument("--input-corrections", help="csv file for host mapping from BIOSAMPLE_Host")
    parser.add_argument("--outfile")

    return parser

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    ARGS = get_parser().parse_args()

    # Define the task
    task_name = "Manually inspection of locations"
    task = f"""This is the task: Manually inspection of the locations:
        1) {ARGS.input_corrections}
         1.1) For each 'input_location' find the appropiate name or coordinates and add it to 
            'output_locations' column.
            - Whenever possible, name is preferred.
            - Coordinates must follow the format lat;lng and the string 'coordinates' 
                must be added to 'notes' column.
            - Name format: 'Country,State,City'

        2) {ARGS.input_locations}
         2.1) Select the new entries of the mapping (filter by version)
         2.2) Select the entries of type == 'name' (we assume coordinates are correct)
         2.3) Evaluate if the output_location is correct. To simplify the process, filter by 
            dist_km (> 20) AND similarity_locs (< 80). We will filter places not sharing the
            same name and far away from each other. These parameters have been calculated
            comparing Google vs Nominatim (openstreet) results.

            - If the location is correct, do nothing
            - If the location is NOT correct, search for the appropiate name or coordinates 
                and add and entry to {ARGS.input_corrections}. 
                    
        NOTE:
            [1] Check the corrected names in https://nominatim.openstreetmap.org/ui/search.html      
            [2] You can find some useful functions for the manual curation process in 
                utils_locations:
                    - 'compute_middle_point'
                    - 'reverse_output_location'
                """

    # Ask the user if the task has been completed
    print(task)
    user_input = input(f"Have you completed the {task_name} ? (yes/no): ").strip().lower()

    # Check the user's response
    if user_input == "yes":
        print(f"Great job! You've completed '{task_name}'.")
        
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