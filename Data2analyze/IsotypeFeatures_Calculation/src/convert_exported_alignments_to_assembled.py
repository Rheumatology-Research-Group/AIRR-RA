#!/usr/bin/env python3

import re
import csv
import sys
import click
import pandas as pd
from pathlib import Path
from loguru import logger
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA


IMGT_REGIONS = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']


class ExportedAlignmentsSample:
    """ File format from mixcr's exportAlignments """

    def __init__(self, file: str):
        self.file = file
        self.df = pd.read_csv(file, delimiter='\t', low_memory=False, quoting=csv.QUOTE_ALL)
        imputed_cols = {col: re.sub('imputed', '', col, flags=re.IGNORECASE)
                        for col in self.df.columns if re.findall('imputed', col, flags=re.IGNORECASE)}
        self.df.rename(columns=imputed_cols, inplace=True)
        self.sample = Path(file).stem
    

    def as_assembled(
        self,
        region_start: str,
        region_end: str,
        species: str,
        is_igh=False,
        keep_non_productive=False
    ) -> pd.DataFrame:
        assert species in ['human', 'mouse']

        df = self.df.copy()
        if region_start == 'full-seq' and region_end == 'full-seq':
            df['aaSeqFR4'] = df['aaSeqFR4'].str.replace('_$', '', regex=True)
            df['v_read_begin'] = df['allVAlignments'].apply(lambda x: int(x.split(';')[0].split('|')[3]))
            df['full_nuc'] = df.apply(lambda row: row['targetSequences'][row['v_read_begin']:], axis=1)
            df['full_pep'] = ""
        elif region_start == 'full-seq' or region_end == 'full-seq':
            sys.exit('If using full-seq, must set both region-start and region-end to full-seq')
        else:
            if region_start not in IMGT_REGIONS:
                sys.exit(f'Invalid sequencing start. Must be in {IMGT_REGIONS}')
            if region_end not in IMGT_REGIONS:
                sys.exit(f'Invalid sequencing end. Must be in {IMGT_REGIONS}')

            regions = IMGT_REGIONS[IMGT_REGIONS.index(region_start):IMGT_REGIONS.index(region_end) + 1]
            required_cols = [f'aaSeq{reg}' for reg in regions] + [f'nSeq{reg}' for reg in regions]
            len1 = len(df)
            df = df[df[required_cols].notnull().all(1)]
            len2 = len(df)
            logger.info(f'{self.sample}: {(len2 / len1) * 100.0:.6g}% of sequences were complete.')

            df['aaSeqFR4'] = df['aaSeqFR4'].str.replace('_$', '', regex=True)
            df['full_pep'] = df[[f'aaSeq{reg}' for reg in regions]].apply(''.join, axis=1)
            df['full_nuc'] = df[[f'nSeq{reg}' for reg in regions]].apply(''.join, axis=1)
            if not keep_non_productive:
                # Filter cdr3s with stop codons
                df = df[~df['full_pep'].str.match('.*[\*_]')]
                len3 = len(df)
                logger.info(f'{self.sample}: {(len3 / len1) * 100.0:.6g}% of complete sequences were productive.')

        df.rename(columns={'aaSeqCDR3': 'cdr3_pep', 'nSeqCDR3': 'cdr3_nuc'}, inplace=True)

        df['V'] = df['allVHitsWithScore'].apply(pick_best_vj_allele)
        df['J'] = df['allJHitsWithScore'].apply(pick_best_vj_allele)

        if is_igh:
            df['C'] = df['allCHitsWithScore'].apply(lambda x: pick_best_isotype(x, species))
            if species == 'human':
                df.loc[df['C'] == 'IGHE', 'C'] = df[df['C'] == 'IGHE']['targetSequences'].apply(lambda x: mixcr_ige_hotfix(x, species))
        else:
            df['C'] = df['allCHitsWithScore'].astype(str).apply(lambda x: re.sub('\*\d+.*$', '', x.split(',')[0]))

        df['C'].fillna('-', inplace=True)

        df = df[['V', 'J', 'C', 'cdr3_pep', 'cdr3_nuc', 'full_pep', 'full_nuc']].copy()
        df['freq'] = 1
        df = df.groupby(['V', 'J', 'C', 'cdr3_pep', 'cdr3_nuc', 'full_pep', 'full_nuc']).sum()
        df.sort_values('freq', ascending=False, inplace=True)

        return df


def pick_best_vj_allele(x):
    if pd.isnull(x):
        return x
    return re.sub('\*\d+\([\d.]+\)', '', x.split(',')[0])


HUMAN_C_REF = {
    'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGC': 'IGHG12',
    'CACCCACCAAGGCTCCGGATGTGTTCCCCATCAT': 'IGHD',
    'GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGT': 'IGHM',
    'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAG': 'IGHA',
    'CCTCCACACAGAGCCCATCCGTCTTCCCCTTGAC': 'IGHE',
    'CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGC': 'IGHG34'
}

# Our primers cannot differentiate all subisotypes of IgG,
# only the difference between IgG1/2 and IgG3/4, so we need
# to combine them. Also, IGHGP contains the same primer sequence
# as IgG1/2, so it is most likely that those calls will actually
# be IgG1/2.
MIXCR_HUMAN_ISOTYPE_COMBOS = {
    'IGHA1-IGHA2': 'IGHA',
    'IGHG1-IGHG2-IGHGP': 'IGHG12',
    'IGHG3-IGHG4': 'IGHG34',
    'IGHG1-IGHG2': 'IGHG12',
    'IGHGP': 'IGHG12',
    'IGHA1': 'IGHA',
    'IGHA2': 'IGHA',
    'IGHG1': 'IGHG12',
    'IGHG2': 'IGHG12',
    'IGHM': 'IGHM',
    'IGHD': 'IGHD',
    'IGHE': 'IGHE'
}


def pick_best_isotype(x, species: str) -> str:
    if pd.isnull(x):
        return x

    hits = [re.findall('(IGH.*)\*.*\((.*)\)', hit)[0] for hit in x.split(',')]
    hits = [(iso, float(score)) for iso, score in hits]
    hits.sort(key=lambda x: x[1], reverse=True)
    hits = [iso for iso, score in hits if score == hits[0][1]]
    hits.sort()

    best = MIXCR_HUMAN_ISOTYPE_COMBOS.get("-".join(hits))
    return best


def mixcr_ige_hotfix(data_seq, species: str) -> str:
    """ We have found a common occurence of C regions that have an IgG1/2 3' end (closest
    to the J region) and IgE 5' end. This is probably a result of an IgE primer picking
    up IgG, which is not unexpected given how common IgG is and rare IgE is.

    This "hotfix" realigns the seqence to a the C reference sequences and picks
    the best isotype.
    """

    # Making sure we don't align to V or J region of sequence
    # Update: I checked and pretty much all C regions are ~20-41 bp, so 50 is fine to be safe
    data_seq = data_seq[-50:]

    alns = [(iso, local_pairwise_align_ssw(DNA(data_seq), DNA(ref_seq))) for ref_seq, iso in HUMAN_C_REF.items()]
    best_iso = max(alns, key=lambda x: x[1][1])[0]
    return best_iso


@click.command()
@click.argument('file', type=click.Path(exists=True))
@click.option('--region-start', '-r', type=click.Choice(['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4', 'full-seq']), required=True)
@click.option('--region-end', '-e', type=click.Choice(['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4', 'full-seq']), required=True, default='FR4')
@click.option('--output', '-o', type=click.Path(), required=False)
@click.option('--is-igh', is_flag=True)
@click.option('--keep-non-productive', is_flag=True)
@click.option('--species', '-s', type=click.Choice(['human', 'mouse']), required=True)
def main(file, region_start, region_end, output, is_igh, keep_non_productive, species):

    sample = ExportedAlignmentsSample(file)
    df = sample.as_assembled(
        region_start,
        region_end,
        species,
        is_igh=is_igh,
        keep_non_productive=keep_non_productive
    )

    if output:
        df.to_csv(output)
    else:
        df.to_csv(sys.stdout)



if __name__ == "__main__":
    main()