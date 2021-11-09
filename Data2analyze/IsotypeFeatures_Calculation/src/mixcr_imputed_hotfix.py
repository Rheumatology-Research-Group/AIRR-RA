#!/usr/bin/env python3

from Bio.Seq import Seq
import click
import pandas as pd
import csv
import re
import sys
import time
import requests
import bs4


VQUEST_URL = 'http://www.imgt.org/IMGT_vquest/analysis'
TAG = 'MYSEQ808'


def get_verified_fr1_pep_from_imgt(species, chain, seq):

    values =  {
        'species': species,
        'receptorOrLocusType': chain,
        'inputType': 'inline',
        'sequences': f'>{TAG}\n{seq}'
    }

    r = requests.post(VQUEST_URL, data=values)
    soup = bs4.BeautifulSoup(r.text, 'html.parser')
    section8 = soup.find(id='sequence1_section8').next_sibling.next_sibling
    lines = section8.getText().split('\n')
    
    for idx, line in enumerate(lines):
        if line.startswith(TAG):
            starting_line = idx
            break
    else:
        raise RuntimeError('Could not find tag in soup')
    
    line1 = lines[starting_line]
    line2 = lines[starting_line + 1]

    part1 = re.split(' {3,}', line1)[1]
    part2 = re.split(' {3,}', line2)[1]

    data_fr1 = part1.split(' ')[0]
    germline_fr1 = part2.split(' ')[0]

    offset = len(re.findall('\.', data_fr1))
    consensus = germline_fr1[:offset].lower().replace('.', '') + re.sub('\.', '', data_fr1)
    return consensus


def translate(seq: str):
    return str(Seq(seq).translate())


@click.command()
@click.argument('file', type=click.Path(exists=True))
@click.option('--reference', '-r', type=click.Path(exists=True), required=True)
@click.option('--output', '-o', type=click.Path(), required=True)
def main(file, reference, output):
    """ Mixcr imputed hotfix

    For paired end, our primers cover part of FR1, but mixcr for
    some reason will only impute FR1 if query_from == 0.

    \b
    reference file generated using:
    repseqio tsv -g VRegion -s human default > reference.tsv
    https://github.com/repseqio/repseqio (library used by mixcr)
    This needs to match the version used to run the data.
    """

    ref = pd.read_csv(reference, delimiter='\t', usecols=['Name', 'Sequence'], index_col='Name')
    ref = ref['Sequence'].to_dict()

    df = pd.read_csv(file, delimiter='\t', low_memory=False, quoting=csv.QUOTE_ALL)

    PARENTHESES_RE = re.compile('\(\d+\)')

    new_rows = []
    for idx, row in df.iterrows():
        # If sequence does not cover CDR1, ignore
        if pd.isnull(row['nSeqImputedCDR1']):
            new_rows.append(row)
            continue

        best_v = PARENTHESES_RE.sub('', row['allVHitsWithScore'].split(',')[0])
        refseq = ref[best_v]
        
        try:
            target_from = int(row['allVAlignments'].split('|')[0])
            query_from = int(row['allVAlignments'].split('|')[3])
            cdr1begin = int(row['refPoints'].split(':')[5])
        except:
            new_rows.append(row)
            continue

        germline_fr1 = refseq[:target_from].lower()
        data_fr1 = row['targetSequences'][query_from:cdr1begin]

        fr1_nuc = None
        fr1_pep = None
        if (len(germline_fr1) + len(data_fr1)) % 3 == 0:
            if len(germline_fr1) % 3 == 0 and len(data_fr1) % 3 == 0:
                fr1_nuc = germline_fr1 + data_fr1
                fr1_pep = translate(germline_fr1).lower() + translate(data_fr1)
                case_no = 0
            else:
                # Allow germline to "steal" n nucleotides from data
                steal_n = 3 - (len(germline_fr1) % 3)
                germline_fr1 = germline_fr1 + data_fr1[:steal_n].lower()
                data_fr1 = data_fr1[steal_n:]
                assert len(germline_fr1) % 3 == 0 and len(data_fr1) % 3 == 0, "error 1"
                fr1_nuc = germline_fr1 + data_fr1
                fr1_pep = translate(germline_fr1).lower() + translate(data_fr1)
                case_no = 1
        else:
            # TODO: This just makes it easier for me
            lower = germline_fr1
            upper = data_fr1

            # TODO: This can be condensed a lot 
            if len(lower) % 3 != 0 and len(upper) % 3 == 0:
                if len(lower) % 3 == 2:
                    # Grab some extra nucleotides from reference to fill in incomplete codon
                    lower = refseq[:target_from + (3 - (len(lower) % 3))].lower()
                    fr1_nuc = lower + upper
                    fr1_pep = translate(lower).lower() + translate(upper)
                    case_no = 2
                elif len(lower) % 3 == 1:
                    lower = lower[:-1]
                    fr1_nuc = lower + upper
                    fr1_pep = translate(lower).lower() + translate(upper)
                    case_no = 3
                else:
                    raise RuntimeError('Error in case 1. Should not happen')
            elif len(lower) % 3 == 0 and len(upper) % 3 != 0:
                
                if len(upper) % 3 == 2:
                    upper = upper[len(upper) % 3:]
                    # Grab some extra nucleotides from reference to fill in incomplete codon
                    lower = refseq[:target_from + 3]
                    fr1_nuc = lower + upper
                    fr1_pep = translate(lower).lower() + translate(upper)
                    case_no = 4
                elif len(upper) % 3 == 1:
                    upper = upper[len(upper) % 3:]
                    fr1_nuc = lower + upper
                    fr1_pep = translate(lower).lower() + translate(upper)
                    case_no = 5
                else:
                    raise RuntimeError('Error in case 2. Should not happen')
            elif len(lower) % 3 != 0 and len(upper) % 3 != 0:
                if len(upper) % 3 == 2:
                    upper = upper[len(upper) % 3:]
                    if len(lower) % 3 == 2:
                        case_no = 6
                        lower = refseq[:target_from + (3 - (len(lower) % 3))]
                    elif len(lower) % 3 == 1:
                        lower = refseq[:target_from + 3 + (3 - (len(lower) % 3))]
                        case_no = 7
                    else:
                        raise RuntimeError('Error in case 3. Should not happen')
                    
                    fr1_nuc = lower + upper
                    fr1_pep = translate(lower).lower() + translate(upper)

                elif len(upper) % 3 == 1:
                    upper = upper[len(upper) % 3:]
                    # Grab some extra nucleotides from reference to fill in incomplete codon
                    if len(lower) % 3 == 2:
                        lower = refseq[:target_from + 3 + (3 - (len(lower) % 3))]
                        case_no = 8
                    elif len(lower) % 3 == 1:
                        lower = refseq[:target_from + (3 - (len(lower) % 3))]
                        case_no = 9
                    else:
                        raise
                    fr1_nuc = lower + upper
                    fr1_pep = translate(lower).lower() + translate(upper)
                else:
                    raise RuntimeError('Error in case 4. Should not happen')
            else:
                raise RuntimeError('Error in case 5. Should not happen')

        row['nSeqImputedFR1'] = fr1_nuc
        row['aaSeqImputedFR1'] = fr1_pep
        new_rows.append(row)

    df2 = pd.DataFrame(new_rows)
    df2.to_csv(output, index=False, sep='\t')


if __name__ == "__main__":
    main()
