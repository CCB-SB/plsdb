
##################################################
# IMPORTS
##################################################
from datetime import datetime

##################################################
# MAIN
##################################################
if __name__ == "__main__":

    # Define the task
    task_name = f"Manually inspection of {snakemake.input.mapping} file"
    task = f"""This is the task: Manually inspection of the mapping file:
        - {snakemake.input.mapping}
         1) Select the entries ECOSYSTEM_check==True.
         2) Evaluate if the entry should contain more than one TAXONOMY_UID
            - If yes, add 'symbiosis' in the ECOSYSTEM_tags
            - If not,
                - Only keep the correct TAXONOMY_UID
                - If necessary, correct ECOSYSTEM_tags accordingly
                - If convient, correct the ecosystem_mapping term that generates the error
         5) Save the document as {snakemake.params.mapping_checked} and continue the pipeline :D
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