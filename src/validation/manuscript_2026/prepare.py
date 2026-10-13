"""Prepare manuscript cohorts without changing denominators or endpoint definitions."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from common import check_inputs, new_directory, require, write_json

MARKERS = ['CCL8', 'CCL13', 'CCL8_CCL13', 'CCL18']
COMPONENTS = ['CCL8+CCL13+ Myeloid DC', 'Myeloid DC', 'B cells', 'T cells', 'Tregs']

def landmark(time, event):
    """Non-durable=1 for EFS event by month 6; exclude earlier censoring."""
    time = pd.to_numeric(pd.Series(time), errors='coerce')
    event = pd.to_numeric(pd.Series(event), errors='coerce')
    valid = time.ge(0) & np.isfinite(time) & event.isin([0, 1])
    return np.where(valid & event.eq(1) & time.le(6), 1,
                    np.where(valid & (time.gt(6) | (time.eq(6) & event.eq(0))), 0, np.nan))

def median_groups(fraction):
    cutoff = fraction.median()
    high = fraction.gt(cutoff).astype(int)
    return cutoff, high, np.where(high, 'High', 'Low')

def normalize_scores(scores):
    require(scores.Mixture.notna().all() and scores.Mixture.is_unique, 'Duplicate/missing Cologne sample key')
    x = scores[COMPONENTS].apply(pd.to_numeric, errors='raise')
    require(np.isfinite(x).all().all() and x.ge(0).all().all(), 'Invalid absolute component scores')
    total = x.sum(axis=1)
    require(total.gt(0).all(), 'Zero absolute-score denominator')
    require(np.allclose(total, scores['Absolute score (sig.score)'], rtol=0, atol=1e-7), 'Absolute-score sum differs')
    return x.div(total, axis=0)

def unique(frame, key, context):
    require(frame[key].notna().all() and frame[key].is_unique, f'Missing/duplicate {context} key')

def prepare(input_dir, output):
    records = check_inputs(input_dir)
    inp = Path(input_dir)
    survival, response, chemo, exclusions, denominator = [], [], [], [], []

    def retain(frame, cohort, endpoint, is_response=False):
        d = frame.copy()
        d['cohort'], d['endpoint'] = cohort, endpoint
        valid = d.fraction.between(0, 1) & np.isfinite(d.fraction)
        if is_response:
            valid &= d.bad.isin([0, 1])
        else:
            valid &= d.event.isin([0, 1]) & d.time_months.ge(0) & np.isfinite(d.time_months)
        for row in d.loc[~valid].itertuples():
            exclusions.append(dict(cohort=cohort, source_id=row.source_id,
                                   marker=row.marker, reason=getattr(row, 'reason', 'Invalid endpoint/fraction')))
        d = d[valid].copy()
        unique(d, 'source_id', cohort+' sample')
        unique(d, 'patient_id', cohort+' patient')
        d['cutoff'], d['high'], d['group'] = median_groups(d.fraction)
        (response if is_response else survival).append(d)

    for cohort in ['monti', 'rchop', 'chop', 'moffitt', 'zuma7']:
        f = pd.read_csv(inp/f'kelvin_{cohort}_fractions.tsv', sep='\t')
        unique(f, 'Mixture', cohort)
        parts = ['CCL8', 'B', 'Myeloid', 'T'] + (['Tregs'] if 'Tregs' in f else [])
        x = f[parts]
        require(np.isfinite(x).all().all() and x.ge(0).all().all(), 'Invalid deconvolution components')
        require(np.allclose(x.sum(axis=1), 1, rtol=0, atol=1e-6), 'Expected relative fractions summing to one')
        f['fraction'] = f.CCL8
        denominator.append(dict(cohort=cohort, definition='All retained reference components', normalized=False))
        if cohort in ['monti', 'rchop', 'chop']:
            c = pd.read_csv(inp/f'{cohort}_clinical.tsv', sep='\t')
            unique(c, 'Array', cohort+' clinical')
            # Historical Figure 5e computes type-7 quartiles BEFORE missing OS omission.
            d = f[['Mixture', 'fraction']].merge(c, left_on='Mixture', right_on='Array', validate='one_to_one')
            d['cohort'], d['source_id'], d['patient_id'] = cohort, d.Mixture, d.Mixture
            require(d.OS_Status.dropna().isin([0, 1]).all(), 'Invalid chemotherapy event coding')
            require(d.OS_Time.dropna().ge(0).all(), 'Negative chemotherapy follow-up')
            chemo.append(d[['cohort', 'source_id', 'patient_id', 'OS_Time', 'OS_Status', 'fraction']])
        else:
            c = pd.read_csv(inp/f'{cohort}_clinical.csv')
            unique(c, 'geo_accession', cohort+' clinical')
            c['source_id'] = c.geo_accession
            c['patient_id'] = c['patient_id'] if 'patient_id' in c else c.geo_accession
            d = c.merge(f[['Mixture', 'fraction']], left_on='source_id', right_on='Mixture', how='left', validate='one_to_one')
            d['marker'], d['pooled'] = 'CCL8_CCL13', True
            if cohort == 'moffitt':
                d['bad'] = d.response_group.map({'DR':0, 'NDR':1})
                d['reason'] = 'Missing/unknown study response'
                endpoint = 'Study DR/NDR (minimum six-month follow-up)'
            else:
                require(d.treatment_arm.eq('Axicabtagene Ciloleucel').all() and d.visit.eq('Screening').all(), 'Unexpected ZUMA-7 arm/visit')
                d['time_months'], d['event'] = d.event_free_survival_months, d.event_free_survival_event
                d['bad'] = landmark(d.time_months, d.event)
                d['reason'] = 'Early censoring, missing EFS or missing fraction'
                endpoint = 'Six-month EFS-derived durability proxy'
            retain(d, cohort, endpoint, True)

    scores = pd.read_csv(inp/'AND_oversampled_scores.tsv', sep='\t')
    shares = normalize_scores(scores)
    scores.to_csv(output/'cologne_absolute_scores_LOCAL_ONLY.csv', index=False)
    require(len(scores) == 24, 'Expected 24 Cologne samples')
    c = pd.read_csv(inp/'stahl_clinical.csv')
    unique(c, 'sample_id', 'Cologne clinical')
    require(set(c.sample_id) == set(scores.Mixture), 'Cologne clinical/sample join is incomplete')
    f = pd.DataFrame({'source_id':scores.Mixture, 'fraction':shares[COMPONENTS[0]]})
    d = c.merge(f, left_on='sample_id', right_on='source_id', validate='one_to_one')
    d['patient_id'], d['marker'], d['pooled'] = d.source_id, 'CCL8_CCL13', True
    require(d.durable_response.eq('DR*').sum() == 2, 'Expected two radiotherapy-associated DR* exclusions')
    d['bad'] = d.durable_response.map({'DR':0, 'NDR':1})
    d['reason'] = np.where(d.durable_response.eq('DR*'), 'DR*: complete remission after additional radiotherapy; excluded', 'Unknown response')
    retain(d, 'stahl', 'Native complete remission at six months', True)
    shares.insert(0, 'source_id', scores.Mixture)
    shares.to_csv(output/'cologne_normalized_shares.csv', index=False)
    denominator.append(dict(cohort='stahl', definition='Sum of all five retained absolute component scores', normalized=True))

    annotations = pd.read_csv(inp/'li_cell_annotations.csv.gz', usecols=['sample_id', 'celltype'], dtype={'sample_id':str})
    require(len(annotations) == 1886312 and annotations.sample_id.nunique() == 232, 'Unexpected LymphoMAP atlas size')
    total = annotations.groupby('sample_id').size().rename('n_total')
    myeloid = annotations[annotations.celltype.isin(['Mac', 'Mo', 'DC', 'pDC'])].groupby('sample_id').size().reindex(total.index, fill_value=0).rename('n_myeloid')
    require(myeloid.sum() == 155180 and myeloid.eq(0).sum() == 1, 'Unexpected myeloid annotations')
    counts = pd.concat([total, myeloid], axis=1).reset_index()
    pd.testing.assert_frame_equal(counts, pd.read_csv(inp/'li_total_cell_counts.csv'), check_dtype=False)
    gates = pd.read_csv(inp/'li_myeloid_gates.csv.gz', dtype={'sample_id':str})
    require(len(gates) == 155180, 'Unexpected marker-gate count')
    for marker in MARKERS:
        require(gates['gate_'+marker].isin([True, False, 0, 1]).all(), 'Invalid marker gate')
    require((gates.gate_CCL8_CCL13 == (gates.gate_CCL8.astype(bool) & gates.gate_CCL13.astype(bool))).all(), 'AND gates inconsistent')
    aggregate = gates.groupby('sample_id')[['gate_'+m for m in MARKERS]].sum()
    require(np.array_equal(gates.groupby('sample_id').size(), myeloid.loc[aggregate.index]), 'Myeloid denominator differs')
    old = pd.read_csv(inp/'li_all_biopsy_counts.csv', dtype={'sample_id':str})
    for marker in MARKERS:
        z = old[old.subset.eq(marker)].set_index('sample_id').sort_index()
        require(z.index.equals(aggregate.index) and np.array_equal(z.n_positive, aggregate['gate_'+marker]), 'Marker counts differ from approved counts')
    book = pd.ExcelFile(inp/'Li_Table_S1.xlsx')
    sheet = next(s for s in book.sheet_names if 'S1B' in s)
    c = pd.read_excel(book, sheet_name=sheet, header=1)
    c.columns = [' '.join(str(s).split()) for s in c.columns]
    flag = c['CAR T immediately following biopsy'].eq(1)
    require(flag.equals(c['Treatment line received immediately following biopsy'].fillna('').str.contains('Axi-cel|Liso-cel', case=False)), 'LymphoMAP treatment flags disagree')
    c = c[flag].copy()
    require(len(c) == 22, 'Expected 22 immediately-pre-CAR-T biopsies')
    c['source_id'] = c['Biopsy Unique ID (patientID-biopsy#)'].astype(str)
    c['patient_id'] = c.source_id.str.rsplit('-', n=1).str[0]
    c['time_months'] = pd.to_numeric(c['EFS (years)'])*12
    c['event'] = c['EFS event (0 = no, 1 = yes)']
    c = c[['source_id', 'patient_id', 'time_months', 'event']].merge(counts, left_on='source_id', right_on='sample_id', validate='one_to_one')
    require(len(c) == 22 and c.n_total.gt(0).all(), 'Incomplete LymphoMAP denominator join')
    for marker in MARKERS:
        z = c.merge(aggregate[['gate_'+marker]], left_on='sample_id', right_index=True, how='left', validate='one_to_one')
        require(z.loc[z['gate_'+marker].isna(), 'n_myeloid'].eq(0).all(), 'Missing nonzero-myeloid numerator')
        z['n_positive'] = z['gate_'+marker].fillna(0).astype(int)
        z['fraction'], z['marker'] = z.n_positive/z.n_total, marker
        retain(z, 'li_cart', 'EFS')
        if marker == 'CCL8_CCL13':
            z['bad'], z['pooled'] = landmark(z.time_months, z.event), True
            z['reason'] = 'Early censoring, missing EFS or missing fraction'
            retain(z, 'li_cart', 'Six-month EFS-derived durability proxy', True)
    cols = ['cohort','endpoint','marker','source_id','patient_id','time_months','event','fraction','cutoff','high','group','n_total','n_myeloid','n_positive']
    s, r = pd.concat(survival), pd.concat(response)
    s.reindex(columns=cols).to_csv(output/'survival_analysis_data.csv', index=False)
    r.reindex(columns=cols+['bad','pooled']).to_csv(output/'response_analysis_data.csv', index=False)
    pd.concat(chemo).to_csv(output/'chemo_analysis_data.csv', index=False)
    pd.DataFrame(exclusions).to_csv(output/'exclusions_LOCAL_ONLY.csv', index=False)
    pd.DataFrame(denominator).to_csv(output/'denominator_audit.csv', index=False)
    counts.to_csv(output/'li_counts_LOCAL_ONLY.csv', index=False)
    roster = r.groupby('cohort').size().to_dict()
    require(roster == {'moffitt':35, 'zuma7':134, 'stahl':22, 'li_cart':22}, 'Response roster differs from approved analysis')
    write_json(output/'preparation_validation.json', dict(status='PASS', response_counts=roster, lymphoma_atlas_cells=len(annotations), myeloid_cells=int(myeloid.sum()), cologne_available=24, cologne_DR_star_exclusions=2, files=records))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args=parser.parse_args()
    prepare(args.input_dir, new_directory(args.output_dir))

if __name__ == '__main__':
    main()
