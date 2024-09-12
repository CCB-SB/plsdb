rule mashdb_sketch:
    input: 
        fasta = "../../results/filtering/deduplication/pls_dedup.fasta", # rules.deduplication.output.fasta
    output:
        join(OUTDIR, "mashdb","plsdb_sketch.msh")
    params:
        cmd = "sketch",
        S = 123, # seed
        k = 21, # k-mer size
        s = 1000, #sketch_size
        o = lambda wildcards, output: splitext(output[0])[0],
        extra = ['-i']
    log:
        join(OUTDIR, "mashdb","plsdb_sketch.log")
    benchmark:
        join(OUTDIR, "mashdb","plsdb_sketch.bench")
    threads: workflow.cores
    wrapper:
        "file:///local/plsdb/master/wrappers/mash"

rule mashdb_distS:
    input:
        ref = rules.mashdb_sketch.output[0],
        query = rules.mashdb_sketch.output[0],
    output:
        join(OUTDIR, "mashdb","sim/mash_dist.tsv")
    params:
        cmd="dist",
        d = 0.00123693, # cutoff
    log:
        join(OUTDIR, "mashdb","sim/mash_dist.log")
    benchmark:
        join(OUTDIR, "mashdb","sim/mash_dist.bench")
    threads: workflow.cores
    wrapper:
        "file:///local/plsdb/master/wrappers/mash"

rule mashdb_sim:
    input: rules.mashdb_distS.output[0]
    output: join(OUTDIR, "mashdb", "sim/plsdb_mashdb_sim.tsv")
    run:
        pairs = set()
        with open(input[0], 'r') as ifile, open(output[0], 'w') as ofile:
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

rule mashdb_dist:
    input:
        ref = rules.mashdb_sketch.output[0],
        query = rules.mashdb_sketch.output[0],
    output:
        join(OUTDIR, "mashdb","db/mash_dist.tsv")
    params:
        cmd="dist",
        extra=['-t']
    log:
        join(OUTDIR, "mashdb","db/mash_dist.log")
    benchmark:
        join(OUTDIR, "mashdb","db/mash_dist.bench")
    threads: workflow.cores
    wrapper:
        "file:///local/plsdb/master/wrappers/mash"