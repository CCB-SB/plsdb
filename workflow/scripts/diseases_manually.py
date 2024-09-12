from datetime import datetime

# Define the task
task_name = f"Manually inspection of Disease Mapping. Generated file: {snakemake.params.mapping_checked}"
task = f"""This is the task: Manually inspection of the Disease Mapping:
    - {snakemake.input.mapping}
        1) Select entries to_check==True. They are duplicates, choose the apropiate entry and delete the other.
        2) Select the new entries of the mapping (filter by version)
        3) Select the entries different from the methods "token_set_ratio" and "wratio" (we assume "exact" is correct)
        4) Ignore wratio results >=95 similarity (they are usually right)
        5) For the resting queries, evaluate if the match is correct:
            - if NO
                - modify the ontology_id accordingly. Look in the Disease Ontology file (do_term.csv/online) and/or 
                    Symp ontology (symp_terms.csv). For more than one id, separate by coma without space.
                - if not ontology_id matching, set to -1.
                - set score = -1
                - set method = manual
        6) If ontology is empty try to fill it manually following the steps off 4.NO section
        7) Complete the tags column choosing one or more tags from the following list:
            [disease,symptom,ecosystem,encoded_trait,healthy]
        7) Save the document as {snakemake.params.mapping_checked} and continue the pipeline :D
            """

# Ask the user if the task has been completed
print(task)

user_input = input(f"Have you completed the {task_name} ? (yes/no): ").strip().lower()

# Check the user's response
if user_input == "yes":
    response = f"Great job! Continue with this amazing work :D"
    
    #save mappings
    # Record datetime dd/mm/YY H:M:S
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    # Create outputfile
    with open(snakemake.output[0], "w+") as out:
        out.write(f"[{dt_string}] The following task has been completed :\n")
        out.write(task)

elif user_input == "no":
    response = f"Don't worry, you can still complete '{task_name}'. Come back again after completion"
else:
    response = "Invalid response. Please enter 'yes' or 'no'."

print(response)
