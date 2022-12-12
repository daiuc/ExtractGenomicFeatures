'''
ExtractGenes: 
    Extract genes

Generally each UTR overlaps with an exon, therefore, 
    1. get utrs
    2. get exons
    3. remove exons that overlaps with utrs
    4. further get 5' utr
    5. further get 3' utr

'''



GTF = '/project2/yangili1/cdai/genome_index/hs38/gencode.v38.primary_assembly.annotation.gtf.gz'
GenomeSize = '/project2/yangili1/cdai/genome_index/hs38/STAR-2.7.7a/chrNameLength.txt'


up, dn = '2000', '1000'
rule ExtractGene:
    message: 'Extract protein coding genes'
    input: 
        gtf = GTF,
        gs = GenomeSize
    output: 
        gene = 'resources/gtftk/proteincoding.gene.gencode.v38.bed',
        gene_ext = 'resources/gtftk/proteincoding.gene_l' + up + '_r' + dn + '.gencode.v38.bed'
    params:
        l = up, # upstream
        r = dn  # downstream
    conda: 'gtftk'
    shell: 
        '''
        gtftk select_by_key \
            -i {input.gtf} \
            -k gene_type -v protein_coding | \
            gtftk select_by_key --select-genes \
                --names gene_id,gene_name \
                --bed-format | \
            sortBed -g {input.gs} -i - | uniq > {output.gene}

        slopBed -s -i {output.gene} -g {input.gs} -l {params.l} -r {params.r} | \
            sortBed -g {input.gs} -i - | uniq > {output.gene_ext}
        
        '''


rule ExtractIntergenic:
    message: 'Extract intergenic region'
    input:
        gene = rules.ExtractGene.output.gene,
        gene_ext = rules.ExtractGene.output.gene_ext,
        gs = GenomeSize
    output:
        intergenic = 'resources/gtftk/proteincoding.intergenic.gencode.v38.bed',
        intergenic_sans_ext = 'resources/gtftk/proteincoding.intergenic_sans_ext.gencode.v38.bed'
    shell:
        '''
        complementBed -i {input.gene} -g {input.gs} | \
            sortBed -g {input.gs} -i - | uniq > {output.intergenic}
        
        complementBed -i {input.gene_ext} -g {input.gs} | \
            sortBed -g {input.gs} -i - | uniq > {output.intergenic_sans_ext}
        '''


rule ExtractTranscript:
    message: 'Extract protein coding transcript'
    input: 
        gtf = GTF,
        gs = GenomeSize
    output: 'resources/gtftk/proteincoding.transcript.gencode.v38.bed'
    conda: 'gtftk'
    shell: 
        '''
        gtftk select_by_key \
            -i {input.gtf} \
            -k gene_type -v protein_coding | \
            gtftk select_by_key --select-transcripts \
                --names gene_id,gene_name \
                --bed-format | \
            sortBed -g {input.gs} -i - | \
            uniq > {output}
        '''


rule ExtractTssTes:
    message: 'Extract TSS and TES from gene'
    input: rules.ExtractGene.output.gene
    output: 
        tss = 'resources/gtftk/proteincoding.tss.gencode.v38.bed',
        tes = 'resources/gtftk/proteincoding.tes.gencode.v38.bed'
    shell:
        '''
        awk 'BEGIN {{OFS = "\t"}};
             {{
                if ($6 == "+")
                    A = $2;
                else
                    A = $3 - 1;
                B = A + 1; 
                print $1, A, B, $4, $5, $6
             }}
            ' {input} > {output.tss}
        
        awk 'BEGIN {{OFS = "\t"}};
             {{
                if ($6 == "+")
                    A = $3 - 1;
                else
                    A = $2;
                B = A + 1;
                print $1, A, B, $4, $5, $6
             }}
            ' {input} > {output.tes}
        '''


rule ExtractExon:
    '''NOTE exons extracted here include UTRs
    '''
    message: 'Extract protein coding genes - exons'
    input: 
        gtf = GTF,
        gs = GenomeSize
    output: 'resources/gtftk/proteincoding.exon.gencode.v38.bed'
    conda: 'gtftk'
    shell: 
        '''
        gtftk select_by_key \
                -i {input.gtf} \
                -k gene_type -v protein_coding | \
            gtftk select_by_key --select-exons \
                --names gene_id,gene_name \
                --bed-format | \
            sortBed -g {input.gs} -i - | \
            uniq > {output}
        '''


rule ExtractUTR:
    '''NOTE UTRs overlap with exons
    '''
    message: 'Extract protein coding genes - UTR'
    input: 
        gtf = GTF,
        gs = GenomeSize
    output: 'resources/gtftk/proteincoding.utr.gencode.v38.bed'
    conda: 'gtftk'
    shell: 
        '''
        gtftk select_by_key \
            -i {input.gtf} \
            -k gene_type -v protein_coding | \
            gtftk select_by_key \
                -k feature -v UTR \
                --names gene_id,gene_name \
                --bed-format | \
            sortBed -g {input.gs} -i - | \
            uniq > {output}
        '''



rule ExtractStartStopCodon:
    '''NOTE UTRs overlap with exons
    '''
    message: 'Extract protein coding genes - start/stop codon'
    input: 
        gtf = GTF,
        gs = GenomeSize
    output: 
        start_codon = 'resources/gtftk/proteincoding.start_codon.gencode.v38.bed',
        stop_codon = 'resources/gtftk/proteincoding.stop_codon.gencode.v38.bed'
    conda: 'gtftk'
    shell: 
        '''
        gtftk select_by_key \
            -i {input.gtf} \
            -k gene_type -v protein_coding | \
            gtftk select_by_key \
                -k feature -v start_codon \
                --names gene_id,gene_name \
                --bed-format | \
            sortBed -g {input.gs} -i - | \
            uniq > {output.start_codon}
                
        gtftk select_by_key \
            -i {input.gtf} \
            -k gene_type -v protein_coding | \
            gtftk select_by_key \
                -k feature -v stop_codon \
                --names gene_id,gene_name \
                --bed-format | \
            sortBed -g {input.gs} -i - | \
            uniq > {output.stop_codon}
        '''



rule ExonSansUtr:
    '''Strategy:
        - In Gencode, all UTR overlaps with exons. Confirmed by the following command:
            intersectBed  -wa -v -a resources/gtftk/proteincoding.utr.gencode.v38.bed -b resources/gtftk/proteincoding.exon.gencode.v38.bed
        - for exons that overlap with utr: true exon = exon - utr
    '''
    message: 'substract UTR from exon'
    input: 
        exon = rules.ExtractExon.output,
        utr  = rules.ExtractUTR.output,
        gs = GenomeSize
    output: 
        exon = 'resources/gtftk/proteincoding.exon_sansUtr.gencode.v38.bed'
    shell: 
        '''
        subtractBed -F .5 -s -sorted \
                -g {input.gs} -a {input.exon} -b {input.utr} | \
            uniq > {output}
        
        '''




rule ExtractIntron:
    '''
    Intron = gene substract exon. 
    This is correct because all UTR overlap with exon in gencode annotations. 
    '''
    message: 'Extract intron'
    input:
        transcript = rules.ExtractTranscript.output,
        exon = rules.ExtractExon.output,
        utr  = rules.ExtractUTR.output,
        gs = GenomeSize
    output: 'resources/gtftk/proteincoding.intron.gencode.v38.bed'
    shell:
        '''
        subtractBed -F 0.99 -s -sorted \
                -g {input.gs} -a {input.transcript} -b {input.exon} | \
            sortBed -g {input.gs} -i - | \
            subtractBed -F 0.99 -s -sorted \
                -g {input.gs} -a - -b {input.utr} | \
            sortBed -g {input.gs} -i - | \
            uniq > {output}
        '''


rule ExtractPromoterAndTail:
    message: 'Extract promoter and tail region'
    input:
        gene_ext = 'resources/gtftk/proteincoding.gene_l2000_r1000.gencode.v38.bed',
        tss = 'resources/gtftk/proteincoding.tss.gencode.v38.bed',
        tes = 'resources/gtftk/proteincoding.tes.gencode.v38.bed',
        intron = 'resources/gtftk/proteincoding.intron.gencode.v38.bed',
        exon = 'resources/gtftk/proteincoding.exon.gencode.v38.bed',
        utr = 'resources/gtftk/proteincoding.utr.gencode.v38.bed',
        gs = GenomeSize
    output:
        promoter = 'resources/gtftk/proteincoding.promoter.gencode.v38.bed',
        tail = 'resources/gtftk/proteincoding.tail.gencode.v38.bed'
    params:
        tmp1 = 'tmp1', tmp2 = 'tmp2'
    shell:
        '''
        echo extracting promoters ...
        subtractBed -F 0.5 -s -sorted -g {input.gs} -a {input.gene_ext} -b {input.intron} | \
            sortBed -g {input.gs} -i - | \
            subtractBed -F 0.5 -s -sorted -g {input.gs} -a - -b {input.exon} | \
            sortBed -g {input.gs} -i - | \
            subtractBed -F 0.5 -s -sorted -g {input.gs} -a - -b {input.utr} | \
            sortBed -g {input.gs} -i - > {params.tmp1}
        
        slopBed -i {input.tss} -g {input.gs} -b 100 | sortBed -g {input.gs} -i - > {params.tmp2}
        
        intersectBed -u -a {params.tmp1} -b {params.tmp2} > {output.promoter}
            
        rm {params.tmp1} {params.tmp2}


        echo extracting tails ...
        subtractBed -F 0.5 -s -sorted -g {input.gs} -a {input.gene_ext} -b {input.intron} | \
            sortBed -g {input.gs} -i - | \
            subtractBed -F 0.5 -s -sorted -g {input.gs} -a - -b {input.exon} | \
            sortBed -g {input.gs} -i - | \
            subtractBed -F 0.5 -s -sorted -g {input.gs} -a - -b {input.utr} | \
            sortBed -g {input.gs} -i - > {params.tmp1}
        
        slopBed -i {input.tes} -g {input.gs} -b 100 | sortBed -g {input.gs} -i - > {params.tmp2}
        
        intersectBed -u -a {params.tmp1} -b {params.tmp2} > {output.tail}
            
        rm {params.tmp1} {params.tmp2}

        '''






# rule ExtractUtr_5and3_prime:
#     '''
#     5' UTR = UTR that overlaps with start_codon or tss (extend 5bp each side)
#     3' UTR = UTR that overlaps with stop_codon or tes (extend 5bp each side)
#     '''
#     message: "Extract 5' and 3' UTR"
#     input: 
#         utr = rules.ExtractUTR.output,
#         start_codon = rules.ExtractStartStopCodon.output.start_codon,
#         stop_codon  = rules.ExtractStartStopCodon.output.stop_codon,
#         tss = rules.ExtractTssTes.output.tss,
#         tes = rules.ExtractTssTes.output.tes,
#         gs = GenomeSize
#     output: 
#         five  = 'resources/gtftk/proteincoding.utr5p.gencode.v38.bed',
#         three = 'resources/gtftk/proteincoding.utr3p.gencode.v38.bed'
#     shell:
#         '''
#         sortBed -g {input.gs} -i <(cat {input.tss} {input.start_codon}) | \
#             slopBed -i - -g {input.gs} -b 10 | \
#             intersectBed -s -u -sorted \
#                 -g {input.gs} -a {input.utr} -b - |
#             uniq > {output.five}

#         sortBed -g {input.gs} -i <(cat {input.tes} {input.stop_codon}) | \
#             slopBed -i - -g {input.gs} -b 10 | \
#             intersectBed -s -u -sorted \
#                 -g {input.gs} -a {input.utr} -b - |
#             uniq > {output.three}
#         '''


# rule ExtractIntron:  # another way to get intron, this doesn't look right though
#     message: 'Extract intron'
#     input: GTF
#     output: 'resources/gtftk/proteincoding.intron2.gencode.v38.bed'
#     conda: 'gtftk'
#     shell:
#         '''
#         gtftk intronic -i {input} -o {output} \
#             --names gene_id,gene_name 
#         '''

























