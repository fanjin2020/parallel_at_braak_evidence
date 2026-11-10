"""Reconcile authorized ADNI event tables with the fixed two-year endpoint target."""
from datetime import date
import calendar
from pathlib import Path
import re
import pandas as pd

TABLES = ('STUDYSUM', 'REGISTRY', 'ADVERSE', 'RECADV', 'TREATDIS', 'OTDISC')


def day(value):
    text = str(value).strip()
    if not re.fullmatch(r'\d{4}-\d{2}-\d{2}', text):
        return None
    try:
        return date.fromisoformat(text)
    except ValueError:
        return None


def rid(value):
    text = str(value).strip()
    return str(int(float(text))) if re.fullmatch(r'\d+(?:\.0+)?', text) else None


def death_bounds(value):
    text = str(value).strip().lower()
    exact = day(text)
    if exact:
        return exact, exact
    year = re.fullmatch(r'(\d{4})(?:-xx-xx)?', text)
    month = re.fullmatch(r'(\d{4})-(\d{2})(?:-xx)?', text)
    try:
        if year and 1900 <= int(year[1]) <= 2100:
            y = int(year[1]); return date(y,1,1),date(y,12,31)
        if month and 1900 <= int(month[1]) <= 2100 and 1 <= int(month[2]) <= 12:
            y,m = int(month[1]),int(month[2])
            return date(y,m,1),date(y,m,calendar.monthrange(y,m)[1])
    except ValueError:
        pass
    return None


def norm(value):
    return str(value).strip().casefold()


def yes(value):
    return norm(value) in ('1', '1.0', 'yes', 'true')


def code_in(value, code):
    text = str(value).strip()
    if not re.fullmatch(r'\d+(?:\.0+)?(?:\s*[:;,]\s*\d+(?:\.0+)?)*', text):
        return False
    return code in [int(float(x)) for x in re.split(r'[:;,]', text)]


def registry_alive(row):
    phase = row.get('COLPROT', '')
    kind = norm(row.get('VISTYPE', ''))
    source = norm(row.get('RGSOURCE', ''))
    conducted = norm(row.get('RGCONDCT', ''))
    if conducted in ('no', '0', '0.0'):
        return False
    if kind in ('not done', 'telephone neuropathology follow-up'):
        return False
    if phase in ('ADNI3', 'ADNI4') and kind in ('3', '4', '3.0', '4.0'):
        return False
    if phase in ('ADNIGO', 'ADNI2') and kind in ('3', '3.0'):
        return False
    if source in ('participant', 'both', '1', '3', '1.0', '3.0'):
        return True
    return phase == 'ADNI1' and yes(row.get('RGCONDCT', '')) and not source


def load_evidence(folder, cutoff):
    evidence, inventory, ignored = {}, [], []
    for table in TABLES:
        path = Path(folder)/(table+'.csv')
        if not path.is_file():
            inventory.append(dict(table=table, status='not_available', rows=0))
            continue
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if 'RID' not in frame:
            raise ValueError(table+' has no RID column')
        inventory.append(dict(table=table, status='read', rows=len(frame)))
        for row in frame.to_dict('records'):
            key = rid(row['RID'])
            if key is None:
                ignored.append(dict(table=table, reason='invalid_RID')); continue
            if norm(row.get('HAS_QC_ERROR', '')) in ('1', '1.0', 'has qc error'):
                ignored.append(dict(table=table, reason='source_QC_error')); continue
            record = evidence.setdefault(key, dict(alive=set(), exact=set(), intervals=set(), dead_by=set(),
                                                   death_flag=False, invalid_death_date=False, reasons=set(), sources=set()))
            phase = row.get('COLPROT', '')
            flag, bounded = False, None
            if table == 'REGISTRY':
                when = day(row.get('EXAMDATE', ''))
                if when and when <= cutoff and registry_alive(row):
                    record['alive'].add(when); record['sources'].add(table)
                if row.get('RGREASON') and when and when <= cutoff:
                    record['reasons'].add(table+':'+row['RGREASON'])
                continue
            if table in ('ADVERSE', 'RECADV'):
                flag = yes(row.get('SAEDEATH' if table == 'ADVERSE' else 'AEHDEATH', ''))
                flag = flag or (table == 'ADVERSE' and norm(row.get('AEOUTCOME', '')) in ('fatal', '5'))
                text = row.get('AEHDTHDT', '')
                interval = death_bounds(text)
                if interval and interval[0] <= cutoff:
                    lo,hi = interval[0],min(interval[1],cutoff)
                    record['intervals'].add((lo,hi)); flag = True
                    if lo == hi:
                        record['exact'].add(lo)
                elif text and not interval:
                    record['invalid_death_date'] = True; flag = True
                elif interval and interval[0] > cutoff:
                    continue
            elif table == 'STUDYSUM':
                flag = norm(row.get('SDPRIMARY', '')) in ('death', '1', '1.0')
                bounded = day(row.get('SDDATE', ''))
                if row.get('SDPRIMARY') and bounded and bounded <= cutoff:
                    record['reasons'].add(table+':'+row['SDPRIMARY'])
            elif table == 'TREATDIS':
                reason = row.get('WDREASON', '')
                flag = norm(reason).startswith('death') or (phase == 'ADNI1' and code_in(reason, 2))
                flag = flag or (phase in ('ADNIGO', 'ADNI2') and code_in(reason, 1))
                bounded = day(row.get('EXAMDATE', ''))
                if reason and bounded and bounded <= cutoff:
                    record['reasons'].add(table+':'+phase+':'+reason)
            elif table == 'OTDISC':
                when = day(row.get('OTDATE', ''))
                if row.get('OTREASON') and when and when <= cutoff:
                    record['reasons'].add(table+':'+row['OTREASON'])
            if bounded and bounded > cutoff:
                continue
            if flag:
                record['death_flag'] = True; record['sources'].add(table)
                if bounded:
                    record['dead_by'].add(bounded)
    return evidence, pd.DataFrame(inventory), pd.DataFrame(ignored, columns=['table', 'reason'])


def classify(upper, cutoff, clinical_alive, event):
    alive = set(clinical_alive) | event.get('alive', set())
    alive = {x for x in alive if x <= cutoff}
    exact = {x for x in event.get('exact', set()) if x <= cutoff}
    bound = {x for x in event.get('dead_by', set()) if x <= cutoff}
    intervals = list(event.get('intervals', set())) + [(d,d) for d in exact]
    if intervals:
        lower_death = max(x[0] for x in intervals)
        upper_death = min([x[1] for x in intervals] + list(bound))
        if lower_death > upper_death or any(x > upper_death for x in alive):
            return 'conflicting_evidence'
        if upper_death <= upper:
            return 'documented_death' if lower_death==upper_death else 'documented_death_interval_before_window_end'
        if lower_death > upper:
            return 'documented_alive'
        if any(x >= upper for x in alive):
            return 'documented_alive'
        return 'unknown_death_timing'
    if bound and any(x > min(bound) for x in alive):
        return 'conflicting_evidence'
    if bound and min(bound) <= upper:
        return 'documented_death_by_window_end'
    if any(x >= upper for x in alive):
        return 'documented_alive'
    if event.get('death_flag', False):
        return 'unknown_death_timing'
    return 'unknown'


def snapshot(folder):
    p = Path(folder)/'DATA_DOWNLOADED_DATE.csv'
    if not p.is_file():
        raise ValueError('Export the local package DATA_DOWNLOADED_DATE; snapshot compatibility must be checked')
    data = pd.read_csv(p, dtype=str, keep_default_na=False)
    dates = {day(v) for v in data.to_numpy().ravel() if day(v)}
    if len(dates) != 1:
        raise ValueError('Expected one unambiguous package data snapshot date')
    return next(iter(dates))


def run(base, tables, work, final, plan, seed, diagnostic=False):
    import outcome_mnar as previous
    import longitudinal as core
    import longitudinal_sensitivity as ls
    from common import write_csv, write_json, sha256
    cutoff, source = ls.cutoff_date(base/'data/private/adni/raw')
    source_date = snapshot(tables)
    if source_date != cutoff:
        raise ValueError(f'R-package snapshot {source_date} differs from retained analysis cutoff {cutoff}; align snapshots before combining')
    evidence, inventory, ignored = load_evidence(tables, cutoff)
    source_inventory = Path(tables)/'inventory.csv'
    if source_inventory.is_file():
        version_info = pd.read_csv(source_inventory, dtype=str)
        if {'table','package','version'}.issubset(version_info):
            inventory = inventory.merge(version_info[['table','package','version']],on='table',how='left',validate='one_to_one')
    write_csv(final/'ADNI_table_inventory.csv', inventory)
    dictionary_path = Path(tables)/'DATADIC.csv'
    if dictionary_path.is_file():
        dictionary = pd.read_csv(dictionary_path, dtype=str, keep_default_na=False)
        if {'TBLNAME','FLDNAME','TEXT','CODE'}.issubset(dictionary):
            fields = [c for c in ['PHASE','TBLNAME','FLDNAME','TEXT','CODE','UNITS'] if c in dictionary]
            write_csv(final/'ADNI_event_variable_dictionary.csv',dictionary.loc[dictionary.TBLNAME.isin(TABLES),fields])
    if len(ignored):
        write_csv(final/'ADNI_ignored_record_counts.csv', ignored.value_counts().rename('n').reset_index())
    if not inventory.status.eq('read').any():
        raise ValueError('No authorized ADNI event tables available')
    work.mkdir(parents=True, exist_ok=True)
    previous.prepare(base, work/'original', plan)
    all_rows = pd.read_csv(work/'original/private/eligibility.csv')
    all_rows['rid'] = all_rows.rid.map(rid)
    visits, _ = core.load_visits(base/'data/private/adni/raw/CDR.csv', set(all_rows.rid), 'cdr')
    statuses, comparisons, details, reasons = [], [], [], []
    for row in all_rows.to_dict('records'):
        key = row['rid']; ev = evidence.get(key, {})
        days = [d for d, values in visits.get(key, {}).items() if len(values) == 1]
        status = classify(day(row['window_end']), cutoff, days, ev)
        statuses.append(status)
        comparisons.append(dict(previous=row['survival'], updated=status,
                                administrative_opportunity=bool(row['administrative_opportunity'])))
        details.append(dict(rid=key, status=status, exact_death_dates=';'.join(map(str, sorted(ev.get('exact', [])))),
                            death_date_intervals=';'.join(str(lo)+'/'+str(hi) for lo,hi in sorted(ev.get('intervals', []))),
                            death_recorded_by_dates=';'.join(map(str, sorted(ev.get('dead_by', [])))),
                            sources=';'.join(sorted(ev.get('sources', [])))))
        for reason in ev.get('reasons', []):
            reasons.append(dict(reason=reason, rid=key))
    all_rows['survival'] = statuses
    private = work/'private'; private.mkdir(exist_ok=True)
    write_csv(private/'eligibility.csv', all_rows)
    write_csv(private/'vital_evidence.csv', details)
    transitions = pd.DataFrame(comparisons).value_counts().rename('n').reset_index()
    write_csv(final/'Survival_classification_changes.csv', transitions)
    if reasons:
        counts = pd.DataFrame(reasons).drop_duplicates().groupby('reason').size().rename('n_participants').reset_index()
        write_csv(final/'Attrition_reason_counts.csv', counts)
    eligible = all_rows.administrative_opportunity.astype(bool)
    selected = all_rows.loc[eligible & all_rows.survival.eq('documented_alive')].reset_index(drop=True)
    flow = pd.read_csv(work/'original/MNAR_sample_flow.csv')
    flow = flow.loc[flow.stage.isin(['baseline_A_positive_eligible','complete_baseline_predictors','full_administrative_opportunity'])].to_dict('records')
    flow += [dict(stage=k, n=len(g)) for k,g in all_rows.loc[eligible].groupby('survival')]
    flow += [dict(stage='selected_survivor_target', n=len(selected)),
             dict(stage='observed_two_year_endpoint', n=int(selected.outcome_cdrsb.notna().sum())),
             dict(stage='missing_two_year_endpoint', n=int(selected.outcome_cdrsb.isna().sum()))]
    write_csv(final/'MNAR_sample_flow.csv', flow)
    write_json(final/'ADNI_vital_provenance.json', dict(snapshot=str(source_date), cutoff=str(cutoff),
        downloads=False, inputs=[dict(file=p.name, sha256=sha256(p)) for p in sorted(Path(tables).glob('*.csv'))],
        rules='Exact event death dates retained. Year/month-only death dates are intervals, never imputed exact days. '
              'Disposition dates are upper bounds, never substituted for death dates. '
              'Conflicts and unresolved timing are counted. Participant/both REGISTRY contacts are distinguished from partner-only contacts.',
        estimand=plan, omitted_tables=inventory.loc[inventory.status.ne('read'),'table'].tolist()))
    result = previous.analyze(selected, work/'mnar', plan, seed, diagnostic)
    for p in (work/'mnar').glob('MNAR_*'):
        if p.is_file():
            import shutil
            shutil.copyfile(p, final/p.name)
    return dict(n_target=len(selected), successful_scenarios=int(result.status.eq('estimated').sum()),
                planned_scenarios=len(result), conditional_survivor_estimand=True)
