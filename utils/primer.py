def primer(inf,outf):

    from Bio import SeqIO
    from StringIO import StringIO

    # check for IUPAC and expand if necessary
    IUPAC = {'R': ['A', 'G'], 'W': ['A', 'T'], 'M': ['A', 'C'], 'S': ['C', 'G'],
             'Y': ['C', 'T'], 'K': ['G', 'T'], 'N':['A','T','G','C'], 'H':['A','C','T'],
             'V': ['A','C','G'], 'D': ['A','G','T'],'B': ['C','G','T'], 'A': ['A'], 'T': ['T'], 'G': ['G'], 'C': ['C']}
    t_seqs = list(SeqIO.parse(open(inf, 'r'), 'fasta'))
    primer_versions= {}
    for seq in t_seqs:
        versions = ['']
        for i in str(seq.seq):
            new_pos = []
            for ver in versions:
                alts = IUPAC[i]
                for alt in alts:
                    new_pos.append(ver+alt)
            versions = new_pos
        primer_versions[seq.id] = versions

    all_versions = []

    for primer in primer_versions:
        i = 0
        for version in primer_versions[primer]:
            all_versions.append(SeqIO.read(StringIO(">%s_%i\n%s"%(primer,i,version)),'fasta'))
            i+=1
    ##add reverse complements for all primers
    with open(outf,'wt') as hndl:
        for ver in all_versions:
            hndl.write(ver.format('fasta'))
            ver.id += '_rc'
            ver.seq = ver.seq.reverse_complement()
            hndl.write(ver.format('fasta'))

