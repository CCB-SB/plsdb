
##################################################
# IMPORTS
##################################################
from datetime import datetime

##################################################
# MAIN
##################################################
if __name__ == "__main__":

    # Define the task
    task_name = "Manually inspection of mapping file"
    task = f"""This is the task: Manually inspection of the mapping file:
        - {snakemake.input.mapping}
         1) Select the new entries of the mapping (filter by version)
         2) Select the entries method=='ete4,exact'. Do a fast looking, but they should be correct.
         3) Select entries method=='X' and evaluate if the tags and host are correct:
            - If yes, eval==1
            - If no, eval==0 and write add the correct tags and
                - If host_associated specified, fill the TAXONOMY_UID (NCBItaxonID)
                - Else, TAXONOMY_UID==-1    
         4) Please, improve the assignment (if convinient)
         5) Save the document and continue the pipeline :D
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
        with open(snakemake.output[0], "w+") as out:
            out.write(f"[{dt_string}] The following task has been completed :\n")
            out.write(task)
        
    elif user_input == "no":
        print(f"Don't worry, you can still complete '{task_name}'. Come back again after completion")
    else:
        print("Invalid response. Please enter 'yes' or 'no'.")