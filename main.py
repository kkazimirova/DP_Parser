import re
import os

QUALITY_LIMIT = 50      # threshold for quality attribute in vcf file. Quality must be > QUALITY_LIMIT
CHAIN_LEN_LIMIT = 10    # threshold for nucleotides chain. Len(chain) must be < CHAIN_LEN_LIMIT
EVENT_MAPPING = {'Homozygous Copy Loss': '0', 'CN Loss': '1', 'CN Gain': '3', 'High Copy Gain': '4'}


class EventRecord:
    def __init__(self, tumor_id, chromosome, region_start, region_stop, event):
        self.tumor_id = tumor_id
        self.chromosome = chromosome
        self.region_start = region_start
        self.region_stop = region_stop
        self.event = event


class MutationRecord:
    def __init__(self, line):
        fields = line.split('\t')
        self.sample_id = fields[0]
        self.chromosome = fields[1]
        self.position = fields[2]
        self.reference = fields[3]
        self.alternative = fields[4]
        self.ad_reference = int(fields[6])
        self.ad_alternative = int(fields[7])
        self.event_code = int(fields[10])


"""
Creates a support file containing parsed target attributes from vcf file.
As input it must be vcf file
As output it creates supoort text file
"""
def parse_vfc_file(input_file, output_file, contains_control_sample, control_sample_idx):
    output_file = open(output_file, "w")
    output_file.write("sample\tchrom\tposition\tref\talt\tgt\tdp\tref_counts\tvar_counts\tcn_n\n")

    with open(input_file) as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("##"):
                continue
            elif line.startswith("#"):
                biopsy_names = parse_vfc_header(line)
            else:
                parse_vfc_line(biopsy_names, line, output_file, contains_control_sample, control_sample_idx)

    output_file.close()


"""
Returning names of tumor samples (biopsy samples) from vcf file
"""
def parse_vfc_header(line):
    line = re.sub(r".*FORMAT", "", line).strip()
    biopsy_names = line.split("\t")
    return biopsy_names

"""
"""
def parse_vfc_line(biopsy_names, line, output_file, contains_control_sample, control_sample_idx):
    line = line.split("\t")
    chrom = line[0]
    pos = line[1]
    ref = line[3]
    alt = line[4]
    quality = float(line[5])
    samples = line[9:]

    if quality > QUALITY_LIMIT and len(alt) <= CHAIN_LEN_LIMIT and len(ref) <= CHAIN_LEN_LIMIT:
        if contains_control_sample:
            control_sample_values = parse_sample(samples[control_sample_idx])
            if not control_sample_has_mutations(control_sample_values):
                for i in range(len(samples)):
                    sample = samples[i]
                    sample_values = parse_sample(sample)
                    if tumor_sample_has_mutation(sample_values):
                        save_tumor_sample(biopsy_names[i], sample_values, chrom, pos, ref, alt, output_file)
        else:
            for i in range(len(samples)):
                sample = samples[i]
                sample_values = parse_sample(sample)
                if tumor_sample_has_mutation(sample_values):
                    save_tumor_sample(biopsy_names[i], sample_values, chrom, pos, ref, alt, output_file)


def parse_sample(sample):
    sample = sample.replace('.', '-1')  # missing values replaced by -1
    fields = sample.split(':')
    gt = fields[0]
    ad = fields[1]
    ad_fields = list(map(lambda x: int(x), ad.split(',')))
    dp = fields[2]

    parsed_sample_values = {'gt': gt, 'ad': ad_fields, 'dp': dp}
    return parsed_sample_values


def control_sample_has_mutations(control_sample_values):
    if control_sample_values['gt'] == "0/0":
        return False
    else:
        return True


def tumor_sample_has_mutation(sample_values):
    if sample_values['gt'] == '0/0':
        return False
    else:
        try:
            ad_ref = sample_values['ad'][0]
            ad_alt = sample_values['ad'][1]
            if ad_alt > 0:
                return True
        except IndexError:
            pass
    return False


def save_tumor_sample(sample_name, sample_values, chrom, pos, ref, alt, output_file):
    output_line = ''
    # output_file.write("sample\tchrom\tposition\tref\talt\tgt\tdp\tref_counts\tvar_counts\tcn_n\n")

    for column in (sample_name, chrom, pos, ref, alt, sample_values['gt'], sample_values['dp']):
        output_line += column
        output_line += '\t'

    for i in range(len(sample_values['ad'])):
        if i != 0:
            output_line += '\t'
        output_line += str(sample_values['ad'][i])

    cn_n = '2'

    output_line += '\t'
    output_line += cn_n
    output_line += '\n'

    if chrom not in ('X', 'Y'):
        output_file.write(output_line)


def read_event_file(events_cache, event_files_dir):
    sample_idx = []

    for file_name in os.listdir(event_files_dir):
        sample_id = file_name.replace('.Tumor_events.txt', '')
        sample_idx.append(sample_id)

        with open(os.path.join(event_files_dir, file_name)) as file:
            lines = file.readlines()
            for line in lines:
                if not (line.startswith("#") or line.startswith("Chromosome")):
                    events_cache = parse_event_file_line(line, sample_id, events_cache)
                else:
                    continue

    return events_cache, sample_idx


def parse_event_file_line(line, tumor_id, events_cache_memory):
    line = line.replace('"', '')
    line = re.sub(r'^.*?chr', 'chr', line)

    chromosome, split_line = line.split(':', maxsplit=1)
    start, split_line = split_line.split('-', maxsplit=1)
    stop, split_line = split_line.split('\t', maxsplit=1)
    event, _ = split_line.split('\t', maxsplit=1)

    chromosome = chromosome.replace("chr", "")
    start = int(start.replace(",", ""))
    stop = int(stop.replace(",", ""))

    event_record = EventRecord(tumor_id, chromosome, start, stop, event)
    events_cache_memory.append(event_record)

    return events_cache_memory


def add_events_to_parsed_vcf(parsed_vcf_file_path, output_dir, tumor_idx, events_cache):
    # output_dir = open(output_dir, "w")
    # output_dir.write("sample\tchrom\tposition\tref\talt\tgt\tdp\tref_counts\tvar_counts\tcn_n\tcn_v\tevent\n")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(parsed_vcf_file_path) as file:
        lines = file.readlines()
        for tumor_id in tumor_idx:

            output_file_name = output_dir + '/' + str(tumor_id)
            output_file = open(output_file_name, "w")
            output_file.write("sample\tchrom\tposition\tref\talt\tgt\tdp\tref_counts\tvar_counts\tcn_n\tcn_v\tevent\n")

            for line in lines:
                if line.startswith("sample"):
                    continue
                else:
                    find_mutation_event(line, tumor_id, events_cache, output_file)

            output_file.close()


def find_mutation_event(line, tumor_id, events_cache, output_file):
    line_splits = line.split('\t')
    sample_id = line_splits[0]
    chrom = line_splits[1]
    position = int(line_splits[2])

    check_position(line, sample_id, chrom, position, tumor_id, events_cache, output_file)


def check_position(line, sample_id, chrom, position, tumor_id, events_cache, output_file):
    if sample_id == tumor_id:
        for events_record in events_cache:
            if events_record.tumor_id == sample_id and \
                    events_record.chromosome == chrom and \
                    events_record.region_start <= position < events_record.region_stop and \
                    events_record.event != 'Allelic Imbalance':  # precoooo

                output_line = line.replace('\n', '')
                output_line += '\t'

                try:
                    output_line += EVENT_MAPPING[events_record.event]
                except KeyError:
                    print('Unknown event found - ', events_record.event)
                    continue

                output_line += '\t'
                output_line += events_record.event
                output_line = output_line + '\n'
                output_file.write(output_line)
    else:
        return None


def generate_pyclone_input(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):

        output_file = os.path.join(output_dir, file_name)
        output_file = open(output_file, "w")
        output_file.write("position_id\tmutation_id\tref_counts\tvar_counts\tcn_n\tcn_v\n")

        with open(os.path.join(input_dir, file_name)) as file:
            lines = file.readlines()

            for i, line in enumerate(lines[1:]):
                fields = line.split('\t')
                ref_counts = fields[7]
                var_counts = fields[8]
                cn_n = fields[9]
                cn_v = fields[10]
                position_id = fields[0] + '_' + fields[1] + '_' + fields[2]  # sample, chromosome, position

                output_line = position_id + '\t' + str(i) + '\t' + ref_counts + '\t' + \
                              var_counts + '\t' + cn_n + '\t' + cn_v + '\n'

                output_file.write(output_line)

        output_file.close()


# generate genotypes according to mutation or reference allele and copy numbers
# NEUPRAVENE
def create_genotypes_file(parsed_file):
    output_file = open("./files/OUTPUT/output_genotypes", "w")
    output_file.write("ID\tCHROM\tPOS\tREF\tALT\tAD_REF\tAD_ALT\tEVENT_CODE\tGENOTYPE\n")

    with open(parsed_file) as file:
        lines = file.readlines()
        for mutation_record in lines[1:]:
            mutation = MutationRecord(mutation_record)
            output_file = generate_genotypes(mutation, output_file)

    output_file.close()


def generate_genotypes(mutation, output_file):
    if mutation.ad_reference > mutation.ad_alternative:
        first = mutation.reference
        second = mutation.alternative
    else:
        first = mutation.alternative
        second = mutation.reference

    if mutation.event_code == 1:
        genotype = first
        line = mutation.sample_id + '\t' + mutation.chromosome + '\t' + mutation.position + '\t' + mutation.reference + '\t'
        line += mutation.alternative + '\t' + str(mutation.ad_reference) + '\t' + str(mutation.ad_alternative) + '\t'
        line += str(mutation.event_code) + '\t' + genotype + '\n'
        output_file.write(line)
    else:
        for i in range(mutation.event_code + 1):
            genotype = ''
            for j in range(i):
                genotype += first

            for j in range(i, mutation.event_code):
                genotype += second

            line = mutation.sample_id + '\t' + mutation.chromosome + '\t' + mutation.position + '\t' + mutation.reference + '\t'
            line += mutation.alternative + '\t' + str(mutation.ad_reference) + '\t' + str(
                mutation.ad_alternative) + '\t'
            line += str(mutation.event_code) + '\t' + genotype + '\n'
            output_file.write(line)

    return output_file


if __name__ == '__main__':
    parse_vfc_file(input_file='./Samples/P12/P12.660085.WGS.HF.Source.somaticVariants.vcf',
                   output_file='./Samples/P12/tmp/parsed_vcf',
                   contains_control_sample=True,
                   control_sample_idx=4
                   )

    events_cache_memory = []
    events_cache_memory, sample_names = read_event_file(events_cache=events_cache_memory,
                                                        event_files_dir='./Samples/P12/events')

    print(sample_names)

    add_events_to_parsed_vcf(parsed_vcf_file_path='./Samples/P12/tmp/parsed_vcf',
                             output_dir='./Samples/P12/tmp/vcf_events',
                             tumor_idx=sample_names,
                             events_cache=events_cache_memory)

    generate_pyclone_input(input_dir='./Samples/P12/tmp/vcf_events',
                           output_dir='./Samples/P12/Pyclone_input')
