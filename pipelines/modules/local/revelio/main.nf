process REVELIO {
    tag "$meta.id"
    label 'process_medium'

    // Wave builds container from conda environment (pysam + samtools)
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.revelio.bam"), emit: bam
    tuple val(meta), path("*.revelio.bam.bai"), emit: bai
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Write revelio.py script (embedded to ensure availability in Wave containers)
    cat > revelio_script.py << 'REVELIO_SCRIPT'
#!/usr/bin/env python

'''
Title: revelio.py
Date: 20190522
Author: Adam Nunn
Description:
    This program takes an input BAM file and masks base positions that are
    potentially influenced by bisulfite-seq as observed from alignments between
    query and reference sequences.
'''

import multiprocessing as mp
import argparse, tempfile, resource, os
import pysam
from array import array

def main(BAM,OUT,QUALITY,THREADS=1,TEMP=None,FASTA=None):
    try: genome = build_genome(FASTA)
    except (TypeError):
        genome = None
        pass

    pool = mp.Pool(int(THREADS))
    jobs = []
    rlimit = resource.getrlimit(resource.RLIMIT_NOFILE)
    rlimit = rlimit[0]
    TDIR = tempfile.mkdtemp(dir=TEMP)

    with pysam.AlignmentFile(BAM, "rb") as original:
        header = original.header.to_dict()
        if "CO" in header: del header["CO"]

        for ref in original.get_index_statistics():
            if ref.mapped > 0:
                if genome: reference = genome[ref.contig]
                else: reference = None
                job = pool.apply_async(worker, (BAM, TDIR, header, ref.contig, reference, QUALITY))
                jobs.append(job)

    bam = recursion(jobs,rlimit,pool,TDIR,bool(genome))
    bam = bam[0].get()
    pool.close()
    pool.join()

    if isinstance(bam, str): os.replace(bam, OUT)
    else: os.replace(bam[0], OUT)
    print("\\n----------------")
    print("Success!\\n")

def build_genome(FASTA):
    genome = {}
    with open(FASTA, 'r') as fasta:
        first = True
        for line in fasta:
            line = line.rstrip()
            if line.startswith(">") and (first == True):
                line = line.split(" ")
                chr = line[0][1:]
                sequence = ''
                first = False
            elif line.startswith(">") and (first == False):
                genome[chr] = sequence.upper()
                line = line.split(" ")
                chr = line[0][1:]
                sequence = ''
            else: sequence += line
        genome[chr] = sequence.upper()
    return genome

def get_aligned_pairs_with_sequence(read,genome):
    if genome:
        pairs = read.get_aligned_pairs()
        return pairs, False
    else:
        try: MD = read.get_tag("MD").upper()
        except:
            print('\\n{}\\n'.format("ERROR: Valid FASTA file must be given if MD tags are missing in BAM file"))
            raise SystemExit(1)
        else:
            read.set_tag("MD",MD,"Z")
            pairs = read.get_aligned_pairs(with_seq=True)
            return pairs, MD

def validate_MD(qseq,MD,CT,count):
    if (("T" in qseq) and (not "C" in MD) and CT) or (("A" in qseq) and (not "G" in MD) and not CT):
        count += 1
    return count

def DoubleMasking(reference,qseq,qual,pairs,CT):
    new_seq, new_qual = "", array('B',[])
    for bp in pairs:
        if bp[0] is None: continue
        elif bp[1] is None:
            new_qual.append(qual[bp[0]])
            new_seq += qseq[bp[0]]
        else:
            if reference: refbase = reference[bp[1]].upper()
            else: refbase = bp[2].upper()

            if qseq[bp[0]] == "T" and refbase == "C" and CT:
                new_seq += "C"
                new_qual.append(0)
            elif qseq[bp[0]] == "A" and refbase == "G" and not CT:
                new_seq += "G"
                new_qual.append(0)
            elif (qseq[bp[0]] == "T" and CT) or (qseq[bp[0]] == "A" and not CT):
                new_seq += qseq[bp[0]]
                new_qual.append(0)
            else:
                new_seq += qseq[bp[0]]
                new_qual.append(qual[bp[0]])
    return new_seq, new_qual

def worker(BAM,TDIR,header,ref,reference,QUALITY):
    name = TDIR + "/" + ref + ".bam"
    with pysam.AlignmentFile(BAM, "rb") as original, pysam.AlignmentFile(name, "wb", header=header) as modified:
        rcount, count = 0, 0
        for alignment in original.fetch(ref):
            rcount += 1
            if (alignment.is_unmapped) or (alignment.is_secondary) or (alignment.is_qcfail): continue
            if (alignment.query_alignment_length == 0): continue

            qseq = alignment.query_sequence
            qual = alignment.query_qualities
            if (qseq == None) or (qual == None): continue

            CT = bool((alignment.is_read1 and not alignment.is_reverse) or (alignment.is_read2 and alignment.is_reverse) or (not alignment.is_paired and not alignment.is_reverse))
            pairs, MD = get_aligned_pairs_with_sequence(alignment,reference)

            if MD:
                try: count = validate_MD(qseq,MD,CT,count)
                except (UnboundLocalError): pass

            new = DoubleMasking(reference,qseq,qual,pairs,CT)
            alignment.query_sequence = new[0]

            if QUALITY: alignment.query_qualities = new[1]
            else: alignment.query_qualities = qual

            modified.write(alignment)
    return name, rcount, count

def merger(TDIR,batch):
    bam = tempfile.mktemp(dir=TDIR)
    arguments = ["-fcp",bam]
    arguments = arguments + batch
    pysam.merge(*arguments)
    return bam

def recursion(jobs,rlimit,pool,TDIR,genome):
    if len(jobs) > 1:
        merges = []
        rcount, count = 0, 0

        for i in range(0, len(jobs), int(rlimit/2)):
            batch = []
            for job in jobs[i:i+int(rlimit/2)]:
                j = job.get()
                if isinstance(j, str): batch.append(j)
                else:
                    batch.append(j[0])
                    rcount += j[1]
                    count += j[2]

            merge = pool.apply_async(merger, (TDIR, batch))
            merges.append(merge)

        if (rcount > 0):
            print("Total alignments read from BAM:\\t\\t{}".format(rcount))
            if not genome:
                perc = (count/rcount)*100
                perc = format(perc, '.2f')
                print("Potentially non-standard MD tags counted: {} (~{}% of total alignments)".format(count, perc))

        jobs = recursion(merges,rlimit,pool,TDIR,genome)
    return jobs

parser = argparse.ArgumentParser(description='Revelio bisulfite masking')
parser.add_argument('inbam', metavar='<BAM>', help='Path to input BAM file')
parser.add_argument('outbam', metavar='<OUT>', help='Path to output BAM file')
parser.add_argument('-f','--fasta', metavar='<FASTA>', help='Path to input FASTA file')
parser.add_argument('-t','--temp', metavar='<TEMP>', help='Path to temp directory')
parser.add_argument('-Q','--quality', dest='QUALITY', help='Apply filter to base qualities', default=False, action='store_true')
parser.add_argument('-T','--threads', dest='THREADS', help='Number of threads', type=int, default=1)

args = parser.parse_args()

if __name__ == '__main__':
    main(args.inbam,args.outbam,args.QUALITY,args.THREADS,args.temp,args.fasta)
REVELIO_SCRIPT

    chmod +x revelio_script.py

    # Add MD tags using samtools calmd (required for revelio)
    samtools calmd -b ${bam} ${fasta} > ${prefix}.calmd.bam 2> calmd.log

    # Index the calmd BAM (required for revelio's get_index_statistics)
    samtools index ${prefix}.calmd.bam

    # Run revelio to mask bisulfite conversions
    python revelio_script.py ${prefix}.calmd.bam ${prefix}.revelio.bam

    # Index the output BAM
    samtools index ${prefix}.revelio.bam

    # Clean up intermediate files
    rm -f ${prefix}.calmd.bam ${prefix}.calmd.bam.bai revelio_script.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revelio: 1.1.0
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.revelio.bam
    touch ${prefix}.revelio.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revelio: stub
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
