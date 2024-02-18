import pysam, itertools
import click
import collections

@click.group()
def cli():
    pass

def allbasesSame(s):
    if len(set(s)) < 3:
        return True
    else:    
        return False

def collapse_cmi(string):
    trimers = list()
    cmis = set()
    while len(string) != 0:
        try:
            # majority vote
            if allbasesSame(string[0:3]):
                base = collections.Counter(string[0:3]).most_common(1)[0][0]
                trimers.append(base)
                string = string[3:]
            else:
            # possible combinations
                trimers.append(set(string[0:3]))
                string = string[3:]
        except:
            print("cmi existing indel or wrong")

    for val in itertools.product(*trimers):
        collapse = ''.join(val)
        cmis.add(collapse)
    return cmis

def cmi_greedy(cmis_sets):
    '''cmis_sets: dictionary
                  key is trimer seq and value is the possible monomer seq combination of the key
                  e.g. {"TTTAAAAAACCCTTTAGTCCCTTTAAACCC": {'TAACTACTAC', 'TAACTGCTAC', 'TAACTTCTAC'}, 
                        "TTTAAAAAACCCTTTGGTCCCTTTAAACCC": {'TAACTGCTAC'}}
        Acommoncmi_inset: int
        left_cmis_set: dictionary
                      keep the remain sets and their element after each cover searching;
                      each set represents a trimer and the element in each set is the 
                      possible monomer seq combination of the trimer
        cmi_setID: dictionary
                   key is the uniq monomer seq among remaining sets for evey cover searching
                   value is all trimers that have that uniq monomer seq.
                   e.g. {"TAACTGCTAC": {"TTTAAAAAACCCTTTAGTCCCTTTAAACCC", "TTTAAAAAACCCTTTGGTCCCTTTAAACCC"}}
        cmi_inset: tuple
                   Two elements. The first one represents the unique monomer seq which covers
                   the most sets of the left_cmis_set. The second one is all sets that is covered by
                   that unique monomer seq
                   e.g.: ("TAACTGCTAC", ["TTTAAAAAACCCTTTAGTCCCTTTAAACCC", "TTTAAAAAACCCTTTGGTCCCTTTAAACCC"])
        greedy_cmis: dictionary
                    key is the trimer seq and value is the unique monomer seq that can 
                    represents the trimer selected by set cover method. 
                    For sets don't share a unique monomer with other sets, always
                    pick the first element in those sets as the one can represents the set.
                    e.g. {"TTTAAAAAACCCTTTAGTCCCTTTAAACCC": {'TAACTGCTAC'}, 
                          "TTTAAAAAACCCTTTGGTCCCTTTAAACCC": {'TAACTGCTAC'}}
        min_sets_cover: int
                    For the example here. Its value is 1
    '''
    
    left_cmis_sets = cmis_sets.copy()
    Acommoncmi_inset = len(cmis_sets)
    
    while Acommoncmi_inset > 1:
        cmi_setID = collections.defaultdict(list)
        for setID, cmi_set in left_cmis_sets.items():
            for cmi in cmi_set:
                cmi_setID[cmi].append(setID)
        
        cmi_setID_sorted = dict(sorted(cmi_setID.items(), key=lambda x: len(x[1]), reverse=True))
        cmi_inset = next(iter(cmi_setID_sorted.items()))

        Acommoncmi_inset = len(cmi_inset[1])

        if Acommoncmi_inset > 1:
            for trimer in cmi_inset[1]:
                cmis_sets[trimer] = {cmi_inset[0]}
                del left_cmis_sets[trimer]
        else:
            break
        if len(left_cmis_sets) == 0:
            break
    else: 
        print('done')
    greedy_cmis = {tri:sulocmi.pop() for tri, sulocmi in cmis_sets.items()}    
    min_sets_cover = len(set(greedy_cmis.values()))
    
    return min_sets_cover, greedy_cmis


@cli.command()
@click.option("-i", "--inbam", help="the input bam file")
@click.option("-t", "--tag", help="The tag add to genes or transcripts")
@click.option("-O", "--outlrbam", help="output the bam with corrected CMI")
@click.option("-o", "--output", help="The output file for genes or transcripts-cell count")
def count(inbam, tag, outlrbam, output, sep = "_"):
    """
    count: int
    n_genetag: int
    i: int
    i_steps: int
    total_count: int
    genetrimers: dictionary
                 e.g. {"geneA": {"TTTAAAAAACCCTTTAGTCCCTTTAAACCC", "TTTAAAAAACCCTTTGGTCCCTTTAAACCC"},
                       "geneB": {"TTTGGGTTTCCCTTTAGTCCCTTTAAACCC", "TTTAAATTTCCCTTTGGTCCCTTTAAACCC"}}

    trimer_with_mono: dictionary
                      e.g. {"geneA_TTTAAAAAACCCTTTAGTCCCTTTAAACCC": {'TAACTGCTAC'},
                            "geneA_TTTAAAAAACCCTTTGGTCCCTTTAAACCC": {'TAACTGCTAC'},
                            "geneB_TTTGGGTTTCCCTTTAGTCCCTTTAAACCC": {'TGTCTACTAC'},
                            "geneB_TTTAAATTTCCCTTTGGTCCCTTTAAACCC" {'TATCTGCTAC'}}
    """
    genetrimers = collections.defaultdict(set)
    count = 0
    n_genetag = 0
    outf = open(output, "w")
    outf.write('%s\t%s\n' % ('gene', 'count'))
    with pysam.AlignmentFile(inbam) as bf:
        for i, r in enumerate(bf):
            genetrimers[tag].add(r.qname.split(sep)[-1])
            n_genetag += 1
        print("The total number of input reads is ", i+1)
        print("The total number of input reads with XT tag is ", n_genetag)

        i_steps = 0
        total_count = 0
        trimer_with_mono = collections.defaultdict(list)
    
    for gene, Trimers in genetrimers.items():
        #keys = str(gene) + '_' + str(trimer)
        if not i_steps % 1000:
            print(i_steps)

        if len(Trimers) == 1:
            count = 1
            total_count += count 
            monocmi = collapse_cmi(Trimers).pop()
            keys = str(gene) + '_' + str(Trimers)
            trimer_with_mono[keys].append(monocmi)
        else:
            corrected_cmis = set(tuple(collapse_cmi(tri)) for tri in Trimers)
            if len(corrected_cmis) == 1:
                count = 1
                total_count += count
                for trimer in Trimers:
                    keys = str(gene) + '_' + str(trimer)
                    trimer_with_mono[keys].append(corrected_cmis)
            else:
                trimer_to_combos = {tri: collapse_cmi(tri) for tri in Trimers}
                count, new_trimer_to_combos = cmi_greedy(trimer_to_combos)
                trimer_with_mono = {str(gene)+"_"+str(tri): mono for tri, mono in new_trimer_to_combos.items()}
                total_count  += count
                
        outf.write('%s\t%s\n' % (gene, count))
        i_steps += 1

    lrbam = pysam.AlignmentFile(inbam, 'rb')
    outbam = pysam.AlignmentFile(outlrbam, 'wb', template=lrbam)
    for lr in lrbam:
        #if lr.has_tag(tag):
            #lrgene = lr.get_tag(tag) 
        lrtrimer = lr.qname.split(sep)[1]
        lrkey = str(tag)+'_'+str(lrtrimer)
        #pdb.set_trace()
        assigned_cmi = trimer_with_mono[lrkey]
        lr.set_tag('CM', assigned_cmi, replace=False)
        outbam.write(lr)
    lrbam.close()
    outbam.close()

    return total_count

if __name__ == "__main__":
    cli()
