from Bio.Blast import NCBIWWW
from Bio import SearchIO

ref_database = {"blastn" : "nt", "blastx" : "nr"}
gene_list = ["accD", "matK", "ndhJ", "rpoB", "rpoC1", "ycf5", "rbcL-a"]


def blast(sample_name, loci_name, loci_seq, service):
    if service=="blastn":
        additional_query =""
    elif service=="blastx":
        additional_query = " AND " + loci_name
    loci_seq = ">%s_%s\n%s"%(sample_name,loci_name,loci_seq)
    result_handle = NCBIWWW.qblast(service, ref_database[service], loci_seq, entrez_query="plants[ORGN]" + additional_query)
    with open(sample_name+"_"+service+".xml", 'w') as out_handle:  # Save the result in the temporary file
        out_handle.write(result_handle.read())
    result_handle.close()


def process_blastx(filename, loci_name):
    blastn_qresult = SearchIO.read(filename, 'blast-xml')
    if len(list(blastn_qresult.hits)) == 0:
        print("No match found")
        return None
    blast_hsp = blastn_qresult[0][0]

    with open("Output_" + filename, 'w') as output_handle:
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("Hit summary\n")
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write(str(blast_hsp))
        output_handle.write('\n')
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("Alignment detail\n")
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("%s - %s\n" % (blast_hsp.aln_all[0][0].seq, blast_hsp.aln_all[0][0].id))
        output_handle.write(blast_hsp.aln_annotation['similarity']+'\n')
        output_handle.write("%s - %s\n" % (blast_hsp.aln_all[0][1].seq, blast_hsp.aln_all[0][1].id))

def process_blastn(filename, loci_name):
    blastn_qresult = SearchIO.read(filename, 'blast-xml')
    if len(list(blastn_qresult.hits)) == 0:
        print("No match found")
        return None
    best_hit = blastn_qresult[0]
    best_hit_species = best_hit.description.split(' ')[0]
    best_hit_genus = best_hit.description.split(' ')[1]
    species_same_genus = set([each_hit.description.split(' ')[0] + " " + each_hit.description.split(' ')[1] for each_hit in
                          blastn_qresult.hits if best_hit_genus in each_hit.description])
    blast_hsp = blastn_qresult[0][0]

    with open("Output_" + filename, 'w') as output_handle:
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("Result Summary\n")
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("Sample name : %s\nLoci name : %s\n" % (filename.split('_')[0], loci_name))
        output_handle.write("Estimated species : %s\n" % (best_hit_genus + " " + best_hit_species))
        output_handle.write("Other specie(s) with same genus : ")
        output_handle.write(str(species_same_genus)+'\n')
        if blast_hsp.hit_strand != blast_hsp.query_strand:
            output_handle.write("Reverse complemented? : Yes\n")
        else:
            output_handle.write("Reverse complemented? : No\n")
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("Hit summary\n")
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write(str(blast_hsp))
        output_handle.write('\n')
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("Alignment detail\n")
        output_handle.write("--------------------------------------------------------------------\n")
        output_handle.write("%s - %s\n" % (blast_hsp.aln_all[0][0].seq, blast_hsp.aln_all[0][0].id))
        output_handle.write(blast_hsp.aln_annotation['similarity']+'\n')
        output_handle.write("%s - %s\n" % (blast_hsp.aln_all[0][1].seq, blast_hsp.aln_all[0][1].id))

def single_loci_blast_wrapper(sample_name, loci_name, loci_seq):
    seq = "".join(loci_seq.split('\n\r'))
    blast(sample_name, loci_name, seq, "blastn")
    process_blastn(sample_name + "_blastn.xml", loci_name)
    if loci_name in gene_list:
        blast(sample_name, loci_name, seq, "blastx")
        process_blastx(sample_name + "_blastx.xml", loci_name)

"""Execution"""

