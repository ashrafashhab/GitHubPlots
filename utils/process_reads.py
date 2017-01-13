# Dependancies
from os import remove, path
import warnings
import glob
from itertools import product
from subprocess import Popen, PIPE
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO, SearchIO
from utils import makedir
from primer import primer

def custom_formatwarning(msg, *a):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = custom_formatwarning

        
def pair_fastq_files(data_dir, 
                     r1_pattern, 
                     r2_pattern,
                     pattern_loc = 'mid', 
                     sample_id_func = None):
    
    """
    data_dir: where to write the processed read files
    r1_pattern: the string in the filename indicating that
    this is read 1
    r2_pattern: the string in the filename indicating that
    this is read 2
    pattern_loc: the position of the read indicator in the file
    name
    sample_func_id: a function that received a filename and 
    returns a sample id
    """
    
    # collect incorporated read 2 files
    seen_r2_files = []
    
    # where is the read indicator in the file name?
    file_glob_pattern = '*{pat}*'
    if pattern_loc == 'mid':
        pass
    elif pattern_loc == 'start':
        file_glob_pattern = '{pat}*'
    elif pattern_loc == 'end':
        file_glob_pattern = '*{pat}'
    else:
        raise IOError('patter_loc takes start, mid or end')
    
    # make sure the data dir path ends with /
    data_dir = data_dir.rstrip('/')+'/'
    
    # make a read 1 file glob pattern
    r1_glob = data_dir + file_glob_pattern.format(pat = r1_pattern)
    
    # group read 1, 2 and sample id
    file_pairs = []
    for f in glob.glob(r1_glob):
        f2 = f.replace(r1_pattern, r2_pattern)
        if not path.exists(f2):
            msg = "No r2 file for %s" % f
            warnings.warn(msg)
            continue
        seen_r2_files.append(f2)
        smpl_id = None
        if sample_id_func:
            smpl_id = sample_id_func(f).replace('_','-')
        file_pairs.append([smpl_id, f, f2])
    
    # check if all read 2 files were caught
    for f2 in glob.glob(r1_glob.replace(r1_pattern, r2_pattern)):
        if not f2 in seen_r2_files:
            warnings.warn("%s is orphan and was ignored" % f2)
    
    return file_pairs




metabeatcline = "metaBEAT.py "
metabeatcline+= "--PCR_primer {outdir}primers.txt "
metabeatcline+= "-Q {outdir}QueryMap.txt "
metabeatcline+= "--n_threads {threads} "
metabeatcline+= "--trim_qual {trim_qual} "
metabeatcline+= "--trim_minlength {trim_minlength} "
metabeatcline+= "--merge --{unmerged}  "
metabeatcline+= "-@ {email} -v "

def process_reads(out_dir,
                  threads,
                  qual_vals,
                  length_vals,
                  email,
                  unmerged="merged_only"):
    
    param_sets = product(qual_vals, length_vals)

    out_dir = out_dir.rstrip('/') + '/'
    
    for param_set in param_sets:

        workdir = out_dir+'minqual%i_minlength%i' % param_set
        if path.exists(workdir):
            warnings.warn("skipping %s, exists" % workdir)
            continue
        makedir(workdir)
        cdcline = "cd %s && {mbcline} && cd .. " % workdir
        mbcline = metabeatcline.format(outdir=out_dir, 
                                       threads=threads, 
                                       trim_qual=param_set[0], 
                                       trim_minlength=param_set[1],
                                       unmerged=unmerged,
                                       email=email)
        cline = cdcline.format(mbcline=mbcline)
        p = Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if len(out) > 0:
            with open(workdir + '/log', 'wt') as hndl:
                hndl.write(out)
        if len(err) > 0:
            print('metaBEAT STDERR')
            print (err)
            #print
            pass

def chimera_detection(out_dir, 
                      qual_vals,
                      length_vals,
                      db
                     ):
    
    param_sets = product(qual_vals, length_vals)
    
    out_dir = out_dir.rstrip('/')+'/'
    
    samples = [l.split('\t')[0] for l in open(out_dir+'QueryMap.txt')]
    
    vcline = "vsearch --uchime_ref {q} --db {db} --nonchimeras {nonch} --chimeras {ch}"
    
    for param_set in param_sets:
        workdir = out_dir+'minqual%i_minlength%i/' % param_set
        for smpl in samples:
            smpldir = workdir+smpl+'/'
            q = smpldir + smpl+'_trimmed.fasta'
            nonch = smpldir + smpl+'_trimmed.non-chimera.fasta'
            ch =  smpldir + smpl+'_trimmed.chimera.fasta'
            cline = vcline.format(q=q, db=db, nonch=nonch, ch=ch)
            p = Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            with open(smpldir+'uchim.log','wt') as hndl:
                hndl.write(out)

# parse metaBEAT_read_stats.csv and *-nonchimeras.fasta

def sequence_counts(basepath):

    data = []
    basepath = basepath.rstrip('/')+'/'
    stats_file = basepath + 'metaBEAT_read_stats.csv'
    for l in open(stats_file, 'r').readlines()[1:]:
        smpl, tot, trimtot, trimpe, trimorph, merged, n1,n2,n3,n4, q \
        = l.rstrip().split(',')
        nonchem = 0
        nonchimf = "%s%s/%s_trimmed.non-chimera.fasta" % (basepath, smpl, smpl)
        if path.exists(nonchimf):
            nonchem = len(list(SeqIO.parse(nonchimf,'fasta')))
        
        line_data = [[smpl, 'half total reads', int(tot)/2.0],
                     [smpl, 'half trimmed reads', int(trimtot)/2.0],
                     [smpl, 'trimmed pairs', int(trimpe)/2.0],
                     [smpl, 'trimmed orphans', int(trimorph)],
                     [smpl, 'merged', int(merged)/2.0],
                     [smpl, 'queries', int(q)],
                     [smpl, 'non-chimeric', int(nonchem)]]
        data += line_data 
            
    return data

def compare_seq_counts_among_param_sets(out_dir, 
                                        qual_vals,
                                        length_vals,
                                        ylim=None):
    
    param_sets = product(qual_vals, length_vals)
    
    out_dir = out_dir.rstrip('/')+'/'
    data = []
    
    for param_set in param_sets:
        workdir = out_dir+'minqual%i_minlength%i/' % param_set
        data += [l+['minqual%i_minlength%i' % param_set] for l in sequence_counts(workdir)]
    
    headers = ['smpl', 'data type', 'count', 'param set']
    
    df =  pd.DataFrame(data, columns=headers)
    #print df
    sns.set(style="ticks")
    # Draw a nested boxplot to show bills by day and sex
    sns.boxplot(x="data type", y="count", hue="param set", data=df, palette="PRGn")
    #sns.despine(offset=10, trim=True)
    plt.xticks(rotation=30)
    if ylim:
        plt.ylim(ylim)
        

def query_length_distribution(out_dir,
                              qual_vals,
                              length_vals):
    
    param_sets = product(qual_vals, length_vals)
    
    out_dir = out_dir.rstrip('/')+'/'
    qlengths = {}
    
    for param_set in param_sets:
        workdir = out_dir+'minqual%i_minlength%i/' % param_set
        stats_file = out_dir + 'QueryMap.txt'
        for l in open(stats_file, 'r').readlines():
            smpl = l.partition('\t')[0]
            
            nonchimqf = "%s%s/%s_trimmed.non-chimera.fasta" % (workdir, smpl, smpl) 
            allqf = "%s%s/%s_trimmed.fasta" % (workdir, smpl, smpl) 
            key = workdir.split('/')[-2]+'_'+smpl
            qlengths[key] = {}
            if path.exists(nonchimqf):
                qlengths[key]['nonchimera'] = [len(r) for r in SeqIO.parse(nonchimqf,'fasta')]
            qlengths[key]['all'] = [len(r) for r in SeqIO.parse(allqf,'fasta')]
    
    sns.set(style="white", palette="muted", color_codes=True)

    # Set up the matplotlib figure
    fig_count = len(qlengths)

    row_num = fig_count/2
    if row_num % 2 > 0:
        row_num += 1
    f, axes = plt.subplots(row_num, 2, figsize=(7, 2.5*row_num), sharex=False)
    sns.despine(left=True)

    row = 0
    col = 0
    for q in qlengths:
        
        # Plot a simple histogram with binsize determined automatically
        a = sns.distplot(qlengths[q]['all'],
                         hist=False,
                         color="red",
                         kde_kws={"shade": True},
                         ax=axes[row, col], axlabel=q, label='all')
        
        if 'nonchimera' in qlengths[q]:
            a = sns.distplot(qlengths[q]['nonchimera'],
                             hist=False,
                             color="blue",
                             kde_kws={"shade": True},
                             ax=axes[row, col], label='non chimera')
            
        if col == 1:
            row += 1
            col = 0
        else:
            col = 1
        a.set_xlim((100, 250))    

    plt.setp(axes, yticks=[])
    plt.tight_layout()

def merge_query_fastas(out_dir,
                       qual_val,
                       length_val,
                       nonchim='auto'):

    out_dir = out_dir.rstrip('/') + '/'
    workdir = out_dir + 'minqual%i_minlength%i/' % (qual_val, length_val)

    fastaglob = None
    splitter = '_trimmed'

    if nonchim == False:
        fastaglob = workdir + '*/*_trimmed.fasta'
    elif nonchim == True:
        fastaglob = workdir + "*/*-non-chimeras.fasta"
        splitter = '-non'
        if len(list(glob.glob(fastaglob))) == 0:
            raise RuntimeError('No nonchimera  fasta files')
    elif  nonchim == 'auto':
        fastaglob = workdir + "*/*-non-chimeras.fasta"
        splitter = '-non'
        if len(list(glob.glob(fastaglob))) == 0:
            fastaglob = workdir + '*/*_trimmed.fasta'
            splitter = '_trimmed'
    else:
        raise IOError('nonchim takse Ture, False or "auto"')


    globalid = 0

    with open(workdir + "seqs.fna",'wt') as hndl:
        for f in glob.glob(fastaglob):
            smpl = f.split('/')[-1].split(splitter)[0].replace('_','-')
            records = SeqIO.parse(f,'fasta')
            for r in records:
                r.id = smpl+'_'+str(globalid)
                globalid += 1
                hndl.write(r.format('fasta'))

def make_qiime_mapping_file(out_dir,
                            qual_val,
                            length_val,
                            metadata=None,
                            fields=None):

    out_dir = out_dir.rstrip('/') + '/'
    mapping = open(out_dir + 'Mapping_%i_%i.txt' % (qual_val, length_val),'wt')
    header = '#SampleID\tBarcodeSequence\tLinkerPrimerSequence'
    if metadata and fields:
        for i in fields:
            header += '\t%s' % i
    mapping.write('%s\n' % header)


    for l in open(metadata,'r'):
        smpl, fstq, r1, r2 = l.rstrip().split('\t')
        mapping.write('%s\tAA\tTT\n' % smpl)

    mapping.close()


def check_primer_presence(out_dir,
                          qual_vals,
                          length_vals,
                          evalue=0.01,
                          span=5):

    param_sets = product(qual_vals, length_vals)

    out_dir = out_dir.rstrip('/') + '/'

    iter_results = {}


    for param_set in param_sets:
        workdir = out_dir + 'minqual%i_minlength%i/' % param_set
        query = "%sprimers.txt" % out_dir
        targets = list(glob.glob("%s*/*_trimmed.fasta" % workdir))
        results = {'Sample':[],'Total queries':[], 'Queries with primer': []}
        for target in targets:
            cline = "makeblastdb -in %s -dbtype nucl" % target
            p = Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            cline = "blastn -query %s -db %s -out temp.xml -outfmt 5 -word_size 7 -penalty -3 -reward 1 -dust no" % (query,target)
            p = Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            queries = SearchIO.parse('temp.xml','blast-xml')
            reads_with_problem = set()
            for q in queries:
                for hit in q:
                    hsp = hit[0]
                    if hsp.evalue < evalue and hsp.aln_span > span:
                        reads_with_problem.add(hsp.hit_id)
            results['Sample'].append(target)
            results['Queries with primer'].append(len(list(reads_with_problem)))
            results['Total queries'].append(len(list(SeqIO.parse(target, 'fasta'))))
            if len(list(reads_with_problem)) > 0:
                warnings.warn("Primer(s) detected in %s" % target)
            remove('temp.xml')
        iter_results[workdir] = results
    return iter_results
