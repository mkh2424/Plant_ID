
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez

class parsed_data(object):
    sample_name = ""
    seq_dic = {}

    def __init__(self, sample_name):
        self.sample_name = sample_name

class output_data(object):
    sample_name = ""
    result_per_loci = []

    def __init__(self, sample_name):
        self.sample_name = sample_name

    def sort_based_on_score(self):
        self.result_per_loci = sorted(self.result_per_loci, key=lambda x:x.score,reverse=True)

class blasted_data(object):  # blast.description class에서 inheritance 할 수 있으면 더 fancy하게 만들수 있을텐데.

    loci_name = ""
    title = ""
    score = 0
    expect = 0
    set_specie_of_same_genus = ()

    def __init__(self,loci_name, title, score, expect, set_specie_of_same_genus):
        self.loci_name = loci_name
        self.title = title
        self.score = score
        self.expect = expect
        self.set_specie_of_same_genus = set_specie_of_same_genus


def input_file(filename):
    with open(filename, 'r') as input_data:
        result = parsed_data(input_data.readline().rstrip('\n'))
        for each_line in input_data:
            if each_line.startswith('>'):
                loci_name = each_line.rstrip('\n').lstrip('>')
                result.seq_dic[loci_name] = ""
            else:
                result.seq_dic[loci_name] += each_line.rstrip('\n')
        return result

def blastn_processing(parsed_data):
    output_storage = output_data(parsed_data.sample_name)

    for loci, seq in parsed_data.seq_dic.items():
        print(loci)
        result_handle = NCBIWWW.qblast("blastn", "nt", seq, entrez_query="plants[ORGN]")
        blastn_record = NCBIXML.read(result_handle)
        if len(list(blastn_record.alignments)) == 0:
            print("No match found")
        else:
            best_match = blastn_record.descriptions[0]
            genus_best_match = best_match.title.split(" ")[1]
            species_same_genus = [each_match.title.split(" ")[1] + " " + each_match.title.split(" ")[2] for each_match in blastn_record.descriptions[1:]
                                  if genus_best_match in each_match.title]
            #if Genus is written in shorthand, should use Entrez search per each match
            print(set(species_same_genus))
            output_storage.result_per_loci.append(blasted_data(loci, best_match.title, best_match.score, best_match.e, set(species_same_genus)))

    output_storage.sort_based_on_score()
    return output_storage

def output_file(filename, blasted_final):
    with open(filename, 'w') as output_handle:
        output_handle.write("%s\n" % blasted_final.sample_name)
        for each_loci_best in blasted_final.result_per_loci:
            output_handle.write("%s %s %d %d\n" %
                                (each_loci_best.loci_name, each_loci_best.title, each_loci_best.score,
                                 each_loci_best.expect))
            output_handle.write("Other speicies with same genus : " + str(each_loci_best.set_specie_of_same_genus) + '\n')

        """
        output_handle.write("Sample name : %s\n" % blasted_final.sample_name) 
        for each_loci_best in blasted_final.result_per_loci:
            output_handle.write("Best match for marker %s : %s\n\t-Score : %d\n\t-E-value : %d\n" %
                                (each_loci_best.loci_name, each_loci_best.title, each_loci_best.score, each_loci_best.expect))
        --> report 용 출력 line
        """




def Get_specie_info(blasted_final):
    locus_organsim = {} #dict for organsim information from Genbank
    Entrez.email = "mkh2424@kaist.ac.kr"
    for each_loci in blasted_final.result_per_loci:
        with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=each_loci.title.split('|')[3]) as handle:
            Seq_record = SeqIO.read(handle,"gb")
            locus_organsim[each_loci.loci_name] = Seq_record.annotations["taxonomy"] #genus 까지
            locus_organsim[each_loci.loci_name].append(Seq_record.annotations["source"].split(" ")[2]) #species를 추가
    return locus_organsim

def loci_matrix(input_locus_organism):
    taxo_checksum = {}
    for name1, taxo1 in input_locus_organism.items():
        taxo_checksum[name1] = {}
        for name2, taxo2 in input_locus_organism.items():
            if name1==name2:
                taxo_checksum[name1][name2] = "="
            elif taxo1[-1] == taxo2[-1]:
                taxo_checksum[name1][name2] = "Species"
            elif taxo1[-2] == taxo2[-2]:
                taxo_checksum[name1][name2] = "Genus"
            elif taxo1[-3] == taxo2[-3]:
                taxo_checksum[name1][name2] = "Family"
            else:
                taxo_checksum[name1][name2] = "None"

    return taxo_checksum

def print_matrix(input_taxo_checksum, filename):
    with open(filename,'a') as output_handle:
        output_handle.write("--Matching comparison matrix--\n")
        output_handle.write("{0:>10}".format(""))
        for loci_name in input_taxo_checksum.keys():
            output_handle.write("{0:>10}".format(loci_name))
        output_handle.write('\n')
        for loci1 in input_taxo_checksum.keys():
            output_handle.write("{0:>10}".format(loci1))
            for loci2 in input_taxo_checksum.keys():
                output_handle.write("{0:>10}".format(input_taxo_checksum[loci1][loci2]))
            output_handle.write('\n')

def multiple_loci_blast_wrapper(seq):
    with open("input_test.fasta", 'w') as input_file_maker:
        input_file_maker.write(seq)
    Parse()
    input_parsed_data = input_file('parsed_test.fasta')
    output_blasted_data = blastn_processing(input_parsed_data)
    output_file("output_test.txt", output_blasted_data)
    print_matrix(loci_matrix(Get_specie_info(output_blasted_data)), "output_test.txt")

def Parse():
    seq_temp = ""
    with open("input_test.fasta", 'r') as input_data:
        with open("parsed_test.fasta", 'w') as output_data:
            output_data.write(input_data.readline().rstrip('\n\r')+'\n')
            for line in input_data:
                output_data.write(input_data.readline().rstrip('\n\r') + '\n')



"""

with open("output_blast_temp.txt", 'wb') as handle:
    pickle.dump(output_blasted_data,handle)

with open("output_blast_temp.txt",'rb') as input_handle:
    output_blasted_data = pickle.load(input_handle)

output_file("output_test.txt", output_blasted_data)
print_matrix(loci_matrix(Get_specie_info(output_blasted_data)), "output_test.txt")


print(output_blasted_data.sample_name)
for each in output_blasted_data.result_per_loci:
    print(each.score)

"""
