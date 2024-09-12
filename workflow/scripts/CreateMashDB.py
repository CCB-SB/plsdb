import os

def run_cmd(cmd, split = True, shell=True):
    import shlex

    """
    Run given CMD
    :param cmd: CMD (str)
    :param split: If the command has to be splitted. Set to False if the command has pipes '|' (bool)
    :param shell: If the command requires shell interpreter. Set to True if the command has pipes '|' (bool)
    :return: cmd (str), status (int), stdout + stderr (str)
    """
    import subprocess

    if split:
        cmd = shlex.split(cmd)

    # check = True, raise an exception if process fail
    try:
        p = subprocess.run(cmd, check=True, 
                        shell=shell, encoding="utf-8",
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE)

        p_stderr = p.stderr
        p_stdout = p.stdout
    
    except subprocess.CalledProcessError as exc:
        print(
            f"Process failed because did not return a successful return code. "
            f"Returned {exc.returncode}\n{exc} \n")
        raise exc

    return cmd, p_stdout, p_stderr

#
# HERE
#

os.makedirs(snakemake.output.DIR, exist_ok=True)

# Sketch
cmd=f"""
mash sketch -S 123 -k 21 -s 1000 -p {snakemake.threads} \
    -o {os.path.splitext(snakemake.output.sketch)[0]} \
    -i {snakemake.input[0]}
"""
run_cmd(cmd, split=False, shell=True)
print("Mash sketch correctly created")

# DistS
cmd=f"""
mash dist -d 0.00123693 -p {snakemake.threads} \
    {snakemake.output.sketch} {snakemake.output.sketch} > {snakemake.output.distS}
"""
run_cmd(cmd, split=False, shell=True)
print("Mash distances between closely plasmid correctly calculated.")

# SIM
pairs = set()
with open(snakemake.output.distS, 'r') as ifile, open(snakemake.output.sim, 'w') as ofile:
    for line in ifile:
        sID, qID, dist, pv, sh = line.rstrip('\n').split('\t')
        pair = sorted([sID, qID])
        pair = (pair[0], pair[1])
        # same ID -> skip
        if sID == qID:
            continue
        # group ID already set -> skip
        if pair in pairs:
            continue
        else:
            pairs.add(pair)
            ofile.write(f'{sID}\t{qID}\n')

# Dist
cmd=f"""
mash dist -t -p {snakemake.threads} \
    {snakemake.output.sketch} {snakemake.output.sketch} > {snakemake.output.dist}
"""
run_cmd(cmd, split=False, shell=True)
print("Mash distances correctly calculated.")