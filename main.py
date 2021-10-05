import re
import os

QUALITY_LIMIT = 50
CHAIN_LEN_LIMIT = 10
EVENT_MAPPING = {'Homozygous Copy Loss': '0', 'CN Loss': '1', 'CN Gain': '3', 'High Copy Gain': '4'}


class EventRecord:
    def __init__(self, tumor_id, chromosome, region_start, region_stop, event):
        self.tumor_id = tumor_id
        self.chromosome = chromosome
        self.region_start = region_start
        self.region_stop = region_stop
        self.event = event

class La:
    def __init__(self, line):
        pole = line.split('\t')
        self.sample_id = pole[0]
        # print(self.sample_id)
        self.chromosome = pole[1]
        # print(self.chromosome)
        self.position = pole[2]
        # print(self.position)
        self.reference = pole[3]
        # print(self.reference)
        self.alternative = pole[4]
        # print(self.alternative)
        self.ad_reference = int(pole[6])
        # print(self.ad_reference)
        self.ad_alternative = int(pole[7])
        # print(self.ad_alternative)
        self.event_code = int(pole[10])
        # print(self.event_code)


def parse_vfc_file(file_path):
    output_file = open("./files/OUTPUT/parsed_vcf.txt", "w")
    output_file.write("SAMPLE\tCHROM\tPOS\tREF\tALT\tGT\tAD_REF\tAD_ALT\tDP\n")

    with open(file_path) as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("##"):
                continue
            elif line.startswith("#"):
                list_of_samples, control_sample_positions = parse_vfc_header(line)
            else:
                parse_vfc_line(list_of_samples, control_sample_positions, line, output_file)

    output_file.close()


def parse_vfc_header(line):
    control_sample_positions = []

    line = re.sub(r".*FORMAT", "", line).strip()
    list_of_samples = line.split("\t")

    for i in range(len(list_of_samples)):
        if list_of_samples[i].endswith("C"):
            control_sample_positions.append(i)

    return list_of_samples, control_sample_positions


def parse_vfc_line(list_of_samples, control_sample_positions, line, output_file):
    line = line.split("\t")
    chrom = line[0]
    pos = line[1]
    ref = line[3]
    alt = line[4]
    quality = float(line[5])
    samples = line[9:]

    if quality > QUALITY_LIMIT and len(alt) <= CHAIN_LEN_LIMIT and len(ref) <= CHAIN_LEN_LIMIT:
        for control_sample_id in control_sample_positions:
            control_sample_values = parse_sample(samples[control_sample_id])
            if control_sample_is_without_mutations(control_sample_values):
                for i in range(len(samples)):
                    sample = samples[i]
                    sample_values = parse_sample(sample)
                    if tumor_sample_has_mutation(sample_values):
                        save_tumor_sample(list_of_samples[i], sample_values, chrom, pos, ref, alt, output_file)


def parse_sample(sample):
    sample = sample.replace('.', '-1')  # missing values replaced by -1
    fields = sample.split(':')
    gt = fields[0]
    ad = fields[1]
    ad_fields = list(map(lambda x: int(x), ad.split(',')))
    dp = fields[2]

    parsed_sample_values = {'gt': gt, 'ad': ad_fields, 'dp': dp}
    return parsed_sample_values


def control_sample_is_without_mutations(control_sample_values):
    if control_sample_values['gt'] == "0/0":
        return True
    else:
        return False


def tumor_sample_has_mutation(sample_values):
    if sample_values['gt'] == '0/0':
        return False
    else:
        ad_ref = sample_values['ad'][0]
        ad_alt = sample_values['ad'][1]
        if ad_alt > 0:
            return True
    return False


def save_tumor_sample(sample_name, sample_values, chrom, pos, ref, alt, output_file):
    a = ''

    for i in (sample_name, chrom, pos, ref, alt, sample_values['gt']):
        a += i
        a += '\t'

    for i in range(len(sample_values['ad'])):
        if i != 0:
            a += '\t'
        a += str(sample_values['ad'][i])

    a += '\t'
    a += sample_values['dp']
    a += '\n'

    output_file.write(a)


def read_event_file(events_cache, event_files_dir, tumor_idx):

    for i, file_name in enumerate(os.listdir(event_files_dir)):
        with open(os.path.join(event_files_dir, file_name)) as file:
            lines = file.readlines()
            for line in lines:
                if line.startswith("chr"):
                    events_cache = parse_event_file_line(line, tumor_idx[i], events_cache)
                else:
                    continue

    return events_cache


def parse_event_file_line(line, tumor_id, events_cache_memory):
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


def add_events_to_parsed_vcf(parsed_vcf_file_path, tumor_idx, events_cache_memory):
    output_file = open("./files/OUTPUT/output_with_events.txt", "w")
    output_file.write("SAMPLE\tCHROM\tPOS\tREF\tALT\tGT\tAD_REF\tAD_ALT\tDP\tEVENT\tEVENT_CODE\n")

    with open(parsed_vcf_file_path) as file:
        lines = file.readlines()
        for tumor_id in tumor_idx:
            for line in lines:
                if line.startswith("SAMPLE"):
                    continue
                else:
                    find_mutation_event(line, tumor_id, events_cache_memory, output_file)

    output_file.close()


def find_mutation_event(line, tumor_id, events_chache_memory, output_file):
    line_splits = line.split('\t')
    sample_id = line_splits[0]
    chrom = line_splits[1]
    position = int(line_splits[2])

    check(line, sample_id, chrom, position, tumor_id, events_chache_memory, output_file)


def check(line, sample_id, chrom, position, tumor_id, events_cache_memory, output_file):
    if sample_id == tumor_id:
        a = 1
        for events_record in events_cache_memory:
            if events_record.tumor_id == sample_id:
                if events_record.chromosome == chrom:
                    if events_record.region_start <= position < events_record.region_stop:
                        if events_record.event != 'Allelic Imbalance':
                            line2 = line.replace('\n', '')
                            line2 = line2 + '\t'
                            line2 = line2 + events_record.event
                            line2 = line2 + '\t'
                            line2 = line2 + EVENT_MAPPING[events_record.event]
                            line2 = line2 + '\n'
                            output_file.write(line2)


    else:
        return None


def create_genotypes_file(parsed_file):
    output_file = open("./files/OUTPUT/output_genotypes", "w")
    output_file.write("ID\tCHROM\tPOS\tREF\tALT\tAD_REF\tAD_ALT\tEVENT_CODE\tGENOTYPE\n")

    with open(parsed_file) as file:
        lines = file.readlines()
        for mutation_record in lines[1:]:
            mutation = La(mutation_record)
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
        for i in range(mutation.event_code+1):
            genotype = ''
            for j in range(i):
                genotype += first

            for j in range(i, mutation.event_code):
                genotype += second

            line = mutation.sample_id + '\t' + mutation.chromosome + '\t' + mutation.position + '\t' + mutation.reference + '\t'
            line += mutation.alternative + '\t' + str(mutation.ad_reference) + '\t' + str(mutation.ad_alternative) + '\t'
            line += str(mutation.event_code) + '\t' + genotype + '\n'
            output_file.write(line)

    return output_file


if __name__ == '__main__':
    parse_vfc_file('./files/VFC/P13.WES.Discovery.vcf')

    tumor_idx = ['I062_007.T', 'I062_015.T', 'I062_022.T', 'I062.033.T']

    events_cache_memory = []
    events_cache_memory = read_event_file(events_cache=events_cache_memory,
                                          event_files_dir='./files/EVENTS',
                                          tumor_idx=tumor_idx)

    add_events_to_parsed_vcf(parsed_vcf_file_path="./files/OUTPUT/parsed_vcf.txt",
                             tumor_idx=tumor_idx,
                             events_cache_memory=events_cache_memory)
    create_genotypes_file("./files/OUTPUT/output_with_events.txt")


