# 1st test
#Test. test_type_name = OR Microbiology Culture OR Histopathology 
#OR Test . order_name = 18S RRNA GENE NUCLEIC ACID DETECTION OR  Panfungal Nucleic Acid Test OR Aspergillus NAT OR Aspergillus Nucleic Acid Detection OR ASPPCR OR ASPERGILLUS ANTIGEN (GALACTOMANNAN) OR Galactomannan Serology OR GALA 
#AND Test_Outcome . test_outcome = Positive 
#AND ref_Organism . class_broad = Mould 

import csv
from datetime import datetime,timedelta
import itertools

# NOTE: episode building used to exclude, permanently and for every future
# run, any test whose only current episode link was to a NOT-IFI episode
# (is_ifi_reviewer_agreed = 0). That was removed (see
# Diagnose/Reports/diagnose_patient_regroup_feedback.txt): a NOT-IFI review
# is a judgement about that episode as a whole, not a statement that the
# test result is negative, and a test can legitimately be evidence for more
# than one episode. Episode building now always runs on the full,
# unfiltered test data; review status is preserved separately, per episode,
# when update_episodes.py matches a freshly-built grouping back to its
# existing DB episode by test-id link.

def format_date_dmy(date_str):
    """Convert a date string to 'day-month-year' format if possible."""
    if not date_str:
        return ''
    for fmt in ('%Y-%m-%d', '%d/%m/%Y', '%Y/%m/%d', '%m/%d/%Y', '%d-%m-%Y'):
        try:
            return datetime.strptime(date_str, fmt).strftime('%d-%m-%Y')
        except Exception:
            continue
    return date_str 

def read_csv_to_dicts(filepath):
    data = []
    with open(filepath, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            data.append(row)
    return data

def sort_dicts_by_id_and_date(data, id_field='t1.patient_id', date_field='t1.collection_date', date_format='%Y-%m-%d'):
    def sort_key(item):
        id_value = int(item[id_field])
        if date_format == 'datetime':
            date_value = item[date_field]
        else:
            date_value = datetime.strptime(item[date_field], date_format)
        # Tiebreakers for same patient+date rows — ensures deterministic ordering.
        # 1. Imaging tests (CT=14, PET=15) sort BEFORE non-imaging on the same date.
        #    This preserves the clinical pattern: CT/PET scans establish the
        #    'Unspecified organism' episode anchor first; cultures then link into it
        #    via the Yeast/Mould rules.  Reversing the order causes imaging tests to
        #    fail to find the anchor (it has already been renamed by the culture) and
        #    get orphaned into standalone episodes.
        # 2. test_id then species as final stable tiebreakers within the same type.
        #    t1.test_id may be ';'-joined in later pipeline stages — use 0 as fallback.
        otype    = item.get('t1.order_type_id', '') or ''
        is_imaging = 0 if otype in ('14', '15') else 1   # 0 sorts before 1
        raw_tid  = item.get('t1.test_id', '') or ''
        test_id  = int(raw_tid) if raw_tid.isdigit() else 0
        species  = item.get('t4.species', '') or ''
        return (id_value, date_value, is_imaging, test_id, species)

    return sorted(data, key=sort_key)

def get_epi_org_and_prev(epi_id,epi_org,prev_test,Prev,row,date_row):

    first_date = Prev[prev_test[0]][prev_test[1]]['first_date']
    last_date = Prev[prev_test[0]][prev_test[1]]['date']
    final_infection_site_id = Prev[prev_test[0]][prev_test[1]]['final_infection_site_id']

    epi_org = Prev[prev_test[0]][prev_test[1]]['epi_org']
    row['t1.collection_date'] = last_date.strftime('%Y-%m-%d')
    if 'Unspecified ' in prev_test[1]:
        old_entry = Prev[prev_test[0]][prev_test[1]]
        # SREQ-114: preserve imaging ancestry — if this episode was ever anchored
        # by a CT/PET scan, carry that fact forward so yeast-exclusion guards
        # remain effective even after a non-imaging Unspecified test joins the
        # episode and overwrites order_type_id.
        old_has_imaging = old_entry.get('has_imaging', old_entry.get('order_type_id', '') in ('14', '15'))
        del Prev[prev_test[0]][prev_test[1]]
        Daux = {'date':date_row,'first_date':date_row,'epi_org':epi_org,'class_broad':row['t4.class_broad'],'class_specific':row['t4.class_specific'],'final_infection_site_id':row['final_infection_site_id'],'order_type_id':row['t1.order_type_id'],'has_imaging':old_has_imaging}
        Prev[prev_test[0]][row['t4.species']] = Daux
    elif (prev_test[1].endswith(' species') and
          not row['t4.species'].endswith(' species') and
          'Unspecified' not in row['t4.species']):
        # SREQ-116: upgrade genus-level key (e.g. "Aspergillus species") to the specific
        # name when a more detailed identification arrives.  Without this, a subsequent
        # different specific organism (e.g. "Aspergillus niger" after "Aspergillus fumigatus"
        # merged into an "Aspergillus species" slot) would incorrectly merge via the still-
        # present generic key.  Unlike the Unspecified upgrade, first_date is preserved
        # because the generic test is a real dated result, not a placeholder.
        # Guard: do NOT fire when the incoming species is Unspecified — that would downgrade
        # a genus-level key (e.g. "Aspergillus species") back to "Unspecified organism",
        # exposing the episode as a merge target for unrelated organisms (e.g. Chaetomium).
        old_entry = Prev[prev_test[0]][prev_test[1]]
        del Prev[prev_test[0]][prev_test[1]]
        Daux = {'date':date_row,'first_date':old_entry['first_date'],'epi_org':epi_org,'class_broad':row['t4.class_broad'],'class_specific':row['t4.class_specific'],'final_infection_site_id':row['final_infection_site_id'],'order_type_id':row['t1.order_type_id']}
        Prev[prev_test[0]][row['t4.species']] = Daux

    return(epi_org,epi_id,row,Prev)

def get_epi_and_prev(epi_id,curr_epi,prev_test,Prev,row,date_row):

    day_limit = {'Yeast':90,'Mould':180,'Dermatophyte':180,'Dimorphic fungus':180,'Unknown':180,'Other':180,'':180, 'Unspecified':180}
    first_date = Prev[prev_test[0]][prev_test[1]]['first_date']
    class_broad = Prev[prev_test[0]][prev_test[1]]['class_broad']
    class_specific = Prev[prev_test[0]][prev_test[1]]['class_specific']
    
    if class_specific == 'Cryptococcus and related species':
        limit = 10000
    else:
        limit = day_limit[class_broad]

    if date_row < first_date + timedelta(days=limit):
        ## same episode_org
        curr_epi = Prev[prev_test[0]][prev_test[1]]['epi']
        Prev[prev_test[0]][prev_test[1]]['date'] = date_row
    else:
        curr_epi = epi_id
        epi_id += 1
        Prev[prev_test[0]][prev_test[1]]['date'] = date_row
        Prev[prev_test[0]][prev_test[1]]['first_date'] = date_row
        Prev[prev_test[0]][prev_test[1]]['epi'] = curr_epi
        Prev[prev_test[0]][prev_test[1]]['class_broad'] = row['t4.class_broad']
        Prev[prev_test[0]][prev_test[1]]['order_type_id'] = row['t1.order_type_id']
    return(curr_epi,epi_id,row,Prev)

def find_prev_for_generic(pid,species,Prev):
    sp = species.rsplit(' ',1)
    last_date = None
    sp_out = None
    for psp in Prev[pid]:
        if sp[0] in psp:
            if not last_date or last_date < Prev[pid][psp]['date']:
                last_date = Prev[pid][psp]['date']
                sp_out = psp
        elif species == 'Candida glabrata complex' and psp == 'Nakaseomyces glabratus':
            # Same recency guard as the branch above -- without it, this special
            # case could overwrite a genuinely more recent match found earlier in
            # the loop, with the outcome depending on dict iteration order rather
            # than which candidate is actually most recent.
            if not last_date or last_date < Prev[pid][psp]['date']:
                last_date = Prev[pid][psp]['date']
                sp_out = psp
    if sp_out:
        return (pid,sp_out)
    else:
        return(None)

def find_prev_for_unspecified(pid,species,Prev,current_order_type_id=None):
    last_date = None
    sp_out = None
    for psp in Prev[pid]:
        #DEBUG
        #if pid == '391':
        #    print('Checking unspecified for',pid,psp,species)
        #    print('Prev data:',Prev[pid][psp])
        if not last_date or last_date < Prev[pid][psp]['date']:
            if species == 'Unspecified yeast' and Prev[pid][psp]['class_broad'] != 'Yeast':
                continue
            if species == 'Unspecified mould' and Prev[pid][psp]['class_broad'] != 'Mould':
                continue
            # Prevent imaging tests (PET/CT) from combining with non-Cryptococcus yeast episodes.
            # Direction 1: Current test is imaging, previous is non-Cryptococcus yeast → skip
            if (current_order_type_id in ['14', '15'] and
                Prev[pid][psp]['class_broad'] == 'Yeast' and
                Prev[pid][psp].get('class_specific', '') != 'Cryptococcus and related species'):
                continue
            # Direction 2: Previous test was imaging (or was originally anchored by imaging),
            # current is non-Cryptococcus yeast → skip.
            # Use has_imaging flag (persisted through Unspecified replacements) rather than
            # order_type_id alone, which may have been overwritten by a later non-imaging test.
            _psp_has_imaging = Prev[pid][psp].get('has_imaging', Prev[pid][psp].get('order_type_id', '') in ('14', '15'))
            if (_psp_has_imaging and species == 'Unspecified yeast'):
                continue
            last_date = Prev[pid][psp]['date']
            sp_out = psp
    if sp_out:
        return (pid,sp_out)
    else:
        return(None)

def get_test_order_type_ids(test_file):
    test_order_type_ids = {}
    with open(test_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            test_order_type_ids[row['t1.test_id']] = row['t1.order_type_id']
    return test_order_type_ids

def add_episode_duration(episodes_data, outdir='Data/tmp'):
    """
    Add episode_duration_days column to episodes data.
    Calculates days between first and last test in each episode.
    """
    # Read test dates from tmp_merged2.csv
    test_dates = {}
    with open(f'{outdir}/tmp_merged2.csv', 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            test_id = row['t1.test_id']
            collection_date = row['t1.collection_date']
            if test_id and collection_date:
                test_dates[test_id] = collection_date
    
    # Calculate duration for each episode
    for episode in episodes_data:
        test_ids = episode['t1.test_id'].split(';')
        episode_dates = []
        
        # Collect all dates for tests in this episode
        for test_id in test_ids:
            if test_id in test_dates:
                try:
                    date_obj = datetime.strptime(test_dates[test_id], '%Y-%m-%d')
                    episode_dates.append(date_obj)
                except:
                    continue
        
        # Calculate duration
        if len(episode_dates) > 1:
            episode_dates.sort()
            first_date = episode_dates[0]
            last_date = episode_dates[-1]
            duration_days = (last_date - first_date).days
            episode['episode_duration_days'] = duration_days
        else:
            # Single test or no valid dates
            episode['episode_duration_days'] = 0
    
    return episodes_data

def filter_episodes_nlp(episodes, test_order_type_ids, outdir='Data/tmp'):
    patient_tests = {}
    with open(f'{outdir}/all_db_tests.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            patient_id = row['patient_id']
            test_id = row['test_id']
            order_type_id = row['order_type']
            test_date = datetime.strptime(row['test_date'], '%Y-%m-%d') if row['test_date'] else None
            if patient_id not in patient_tests:
                patient_tests[patient_id] = []
            patient_tests[patient_id].append((test_id, test_date, order_type_id))

    filtered_episodes = []
    for episode in episodes:
        
        diagnosis_date = datetime.strptime(episode['system_diagnosis_date'], '%Y-%m-%d')
        for test_id in episode['t1.test_id'].split(';'):
            order_type_id = test_order_type_ids[test_id]  
            if order_type_id not in ['14', '15']:
                filtered_episodes.append(episode)
                break
        else:
            patient_id = episode['t1.patient_id']
            if patient_id in patient_tests:
                for test_id, test_date, order_type in patient_tests[patient_id]:
                    if order_type not in ['1','2','5','6','8','11']:
                        continue
                    if test_date and -28 <= (diagnosis_date - test_date).days <= 28:
                        filtered_episodes.append(episode)
                        break
    return filtered_episodes

def filter_followup_pet(episodes, test_order_type_ids):
    """
    Merge follow-up imaging episodes into the most recent prior episode from
    the same patient, provided:
      - the current episode date is within 180 days AFTER the prior episode's
        diagnosis date, AND
      - the prior episode is NOT a non-Cryptococcus Yeast episode
        (rule: imaging links for everything except Yeast, excluding Cryptococcus).

    Two tiers of eligibility:

    1. **Imaging-only** episodes (every test is order_type 14 or 15) are
       always eligible for merging into any qualifying prior episode.

    2. **Mixed-Unspecified** episodes (imaging + non-imaging tests where
       class_broad == "Unspecified", i.e. the non-imaging tests did not
       identify an organism) are also eligible.  This handles follow-up
       monitoring that pairs PET/CT scans with histopath, molecular, or
       culture tests that came back negative.

    In both tiers the prior episode must be within 180 days and must not be
    non-Cryptococcus Yeast.
    """
    if not episodes:
        return episodes

    merged_log = []  # collect merge info for reporting

    # Sort episodes by patient_id and system_diagnosis_date
    sorted_episodes = sorted(episodes, key=lambda x: (
        int(x['t1.patient_id']),
        datetime.strptime(x['system_diagnosis_date'], '%Y-%m-%d')
    ))

    filtered_episodes = []

    for current_episode in sorted_episodes:
        patient_id = current_episode['t1.patient_id']
        current_date = datetime.strptime(current_episode['system_diagnosis_date'], '%Y-%m-%d')

        # Get all test IDs for this episode
        test_ids = current_episode['t1.test_id'].split(';')

        # Classify whether the episode contains imaging tests
        has_imaging = any(
            test_order_type_ids.get(tid, '') in ['14', '15']
            for tid in test_ids
        )
        only_imaging = has_imaging and all(
            test_order_type_ids.get(tid, '') in ['14', '15']
            for tid in test_ids
        )

        # Tier 1: imaging-only → always eligible
        # Tier 2: mixed with class_broad "Unspecified" → also eligible
        current_class_broad = current_episode.get('t4.class_broad', '')
        is_mixed_unspecified = (has_imaging and not only_imaging
                                and current_class_broad == 'Unspecified')
        eligible_for_merge = only_imaging or is_mixed_unspecified

        if eligible_for_merge:
            found_previous = False

            # Search backwards through already-processed episodes for this patient
            for prev_episode in reversed(filtered_episodes):
                if prev_episode['t1.patient_id'] != patient_id:
                    continue

                prev_date = datetime.strptime(prev_episode['system_diagnosis_date'], '%Y-%m-%d')
                days_diff = (current_date - prev_date).days

                # Must be within 180 days AFTER the prior episode's diagnosis date
                if not (0 <= days_diff <= 180):
                    continue

                # Skip non-Cryptococcus Yeast episodes (imaging does NOT extend Yeast)
                prev_class_broad    = prev_episode.get('t4.class_broad', '')
                prev_class_specific = prev_episode.get('t4.class_specific', '')
                if (prev_class_broad == 'Yeast' and
                        prev_class_specific != 'Cryptococcus and related species'):
                    continue

                # Tier 2 (mixed-Unspecified) should only merge into a named
                # episode (Mould, Crypto, etc.), not into another Unspecified
                # prior.  This avoids cascading changes in Unspecified-only
                # chains that would displace existing imaging-only merges.
                if not only_imaging and prev_class_broad == 'Unspecified':
                    continue

                # Merge current tests into the prior episode
                prev_test_ids = prev_episode['t1.test_id'].split(';')
                prev_episode['t1.test_id'] = ';'.join(prev_test_ids + test_ids)

                merge_type = 'imaging-only' if only_imaging else 'mixed-unspecified'
                merged_log.append({
                    'patient_id':      patient_id,
                    'imaging_date':    current_date.strftime('%Y-%m-%d'),
                    'days_after_diag': days_diff,
                    'merged_into_diag_date': prev_date.strftime('%Y-%m-%d'),
                    'prev_class_broad': prev_class_broad,
                    'prev_class_specific': prev_class_specific,
                    'merge_type':       merge_type,
                    'imaging_test_ids': ';'.join(test_ids),
                })
                found_previous = True
                break

            if not found_previous:
                filtered_episodes.append(current_episode)
        else:
            # Episode contains non-imaging tests and is not Unspecified — keep as-is
            filtered_episodes.append(current_episode)

    # Print merge report
    if merged_log:
        print(f"\n[filter_followup_pet] Merged {len(merged_log)} episode(s) into prior episodes:")
        print(f"  {'patient_id':<12} {'date':<14} {'days':<8} {'prior_date':<16} {'prior_class':<35} {'type':<20} test_ids")
        print(f"  {'-'*120}")
        for m in merged_log:
            prev_cls = m['prev_class_specific'] if m['prev_class_specific'] != m['prev_class_broad'] else m['prev_class_broad']
            print(f"  {m['patient_id']:<12} {m['imaging_date']:<14} {m['days_after_diag']:<8} {m['merged_into_diag_date']:<16} {prev_cls:<35} {m['merge_type']:<20} {m['imaging_test_ids']}")
    else:
        print("[filter_followup_pet] No imaging episodes merged.")

    return filtered_episodes

def group_episode(outdir='Data/tmp'):
    epi_out = f'{outdir}/epi_table.csv'
    data = read_csv_to_dicts(f'{outdir}/epi_org_table.csv')
    sorted_data = sort_dicts_by_id_and_date(data,date_field='system_diagnosis_date')
    epi_id = 1
    Prev = {}
    out_l = []
    date_format='%Y-%m-%d'
    group_field = 't4.class_broad'
    
    # Get test order type mapping once
    test_order_type_ids = get_test_order_type_ids(f'{outdir}/tmp_merged2.csv')
    for row in sorted_data:
        # Add order_type_id to row by looking it up
        row['t1.order_type_id'] = test_order_type_ids.get(row['t1.test_id'], '')
        
        if row['t4.class_specific'] == 'Cryptococcus and related species':
            group_field = 't4.class_specific'
        else:
            group_field = 't4.class_broad'
        date_row = datetime.strptime(row['system_diagnosis_date'],date_format)
        prev_test = None
        if row['t1.patient_id'] in Prev:
            if row[group_field] in Prev[row['t1.patient_id']]:
                prev_test = (row['t1.patient_id'],row[group_field])
        if prev_test:
            (curr_epi,epi_id,row,Prev) = get_epi_and_prev(epi_id,curr_epi,prev_test,Prev,row,date_row)
        else:
            curr_epi = epi_id
            epi_id += 1
            Daux = {'date':date_row,'first_date':date_row,'epi':curr_epi,'class_broad':row['t4.class_broad'],'class_specific':row['t4.class_specific'],'order_type_id':row['t1.order_type_id']}
            if row['t1.patient_id'] in Prev:
                Prev[row['t1.patient_id']][row[group_field]] = Daux
            else:
                Prev[row['t1.patient_id']] = {row[group_field]:Daux}
        row['episode_id'] = curr_epi
        out_l.append(row)
    # SREQ-115: merge follow-up imaging into prior episodes FIRST, before the
    # SREQ-114 gate drops standalone imaging episodes that lack a nearby supporting
    # test.  Any imaging episode that does NOT match a prior episode within 180 days
    # is left in the list and then evaluated by filter_episodes_nlp below.
    out_l = filter_followup_pet(out_l, test_order_type_ids)
    out_l = filter_episodes_nlp(out_l, test_order_type_ids, outdir)
    
    # Add episode duration before writing to CSV
    out_l = add_episode_duration(out_l, outdir)
    
    with open(epi_out, 'w', newline='', encoding='utf-8') as outfile:
        h = ['t1.patient_id','episode_id','episode_organism_id','t1.test_id','system_diagnosis_date','infection_site_name','t4.species','t4.genus','t4.class_specific','t4.class_broad','system_classification_id','system_ifi_classification','episode_duration_days']
        writer = csv.DictWriter(outfile, fieldnames=h)
        writer.writeheader()
        for r in out_l:
            if 'system_diagnosis_date' in r:
                r['system_diagnosis_date'] = format_date_dmy(r['system_diagnosis_date'])
            # Remove t1.order_type_id before writing (it was only needed for episode grouping)
            if 't1.order_type_id' in r:
                del r['t1.order_type_id']
            writer.writerow(r)

def prev_test_recent(prev_test,Prev,row):
    debug = False
    deb_pat = '391'
    deb_test = '110588'

    date_format='%Y-%m-%d'
    day_limit = {'Yeast':90,'Mould':180,'Dermatophyte':180,'Dimorphic fungus':180,'Unknown':180,'Other':180,'':180, 'Unspecified':180}
    first_date = Prev[prev_test[0]][prev_test[1]]['first_date']
    class_broad = Prev[prev_test[0]][prev_test[1]]['class_broad']
    class_specific = Prev[prev_test[0]][prev_test[1]].get('class_specific', '')
    final_infection_site_id = Prev[prev_test[0]][prev_test[1]]['final_infection_site_id']
    date_row = datetime.strptime(row['t1.collection_date'],date_format)
    prev_order_type = Prev[prev_test[0]][prev_test[1]].get('order_type_id', '')
    
    if (row['t4.species'] == 'Unspecified organism' and 
        (row['t1.order_type_id'] == '14' or row['t1.order_type_id'] == '15')):
        # PET/CT scans with "Unspecified organism" use 28-day window
        daylimit = 28
    elif (prev_test[1] == 'Unspecified organism' and 
          (prev_order_type == '14' or prev_order_type == '15')):
        # Previous test was "Unspecified organism" AND was imaging-only (PET/CT)
        daylimit = 28
    elif (prev_test[1] == 'Unspecified organism' and 
          row['t4.species'] != 'Unspecified organism'):
        # Previous test was "Unspecified organism" but was NOT imaging.
        # Use the INCOMING row's class_broad to determine the window, so that
        # Yeast (90 days) and Mould/Other (180 days) are applied correctly rather
        # than always using 180 days from the Unspecified organism's own class_broad.
        daylimit = day_limit.get(row['t4.class_broad'], 180)
    elif row['t4.class_specific'] == 'Cryptococcus and related species':
        daylimit = 10000
    else:
        # Default: use class_broad time limits (Yeast=90, Mould=180, etc.)
        # This includes: current test is "Unspecified organism" but NOT type 14/15
        daylimit = day_limit[class_broad]

    if debug and row['t1.patient_id'] == deb_pat and row['t1.test_id'] == deb_test:
        print('prev_test:',prev_test)
        print('current:',row['t4.species'])
        print('current type:',row['t1.order_type_id'])
        print('daylimit:',daylimit)
        print('first_date:',first_date)
        print('date_row:',date_row)
        difference = (date_row - first_date).days
        print('difference in days:',difference)
    
    if date_row < first_date + timedelta(days=daylimit):
        return(True)
    else:
        return(False)

def get_prev_test(row,group_field,Prev):

    aux = row[group_field].split(' ',1)
    complex = row[group_field] + ' complex'
    generic = aux[0] + ' species'
    prev_test = None

    debug = False
    deb_pat = '391'
    deb_test = '110588'
    if debug and row['t1.patient_id'] == deb_pat:
        print('Getting prev test for',row['t1.test_id'],row['t1.patient_id'],row[group_field])

    if debug and row['t1.patient_id'] == deb_pat and row['t1.test_id'] == deb_test:
        print('Getting prev test for',row['t1.test_id'],row['t1.patient_id'],row[group_field])
        print('Group field:',row[group_field])
        print('Class broad:',row['t4.class_broad'])
        print('Class specific:',row['t4.class_specific'])
        print('test order type id:',row['t1.order_type_id'])
        if row['t1.patient_id'] in Prev:
            print('Prev:',Prev[deb_pat])
        else:
            print('No previous tests for this patient')

    if row['t1.patient_id'] in Prev:
        if 'Unspecified ' in row[group_field]:
            prev_test = find_prev_for_unspecified(row['t1.patient_id'],row[group_field],Prev,row.get('t1.order_type_id',''))
            if debug and row['t1.test_id'] == deb_test:
                print('Found unspecified',prev_test)
                if prev_test:
                    print('Recent?',prev_test_recent(prev_test,Prev,row))
        if ((not prev_test or not prev_test_recent(prev_test,Prev,row)) and 
            row[group_field] in Prev[row['t1.patient_id']]):
            prev_test = (row['t1.patient_id'],row[group_field])
            if debug and row['t1.test_id'] == deb_test:
                print('Found match',prev_test)
        if ((not prev_test or not prev_test_recent(prev_test,Prev,row)) and
            generic in Prev[row['t1.patient_id']]):
            prev_test = (row['t1.patient_id'],generic)
            if debug and row['t1.test_id'] == deb_test:
                print('Found generic',prev_test)
        if ((not prev_test or not prev_test_recent(prev_test,Prev,row)) and 
            complex in Prev[row['t1.patient_id']]):
            prev_test = (row['t1.patient_id'],complex)
            if debug and row['t1.test_id'] == deb_test:
                print('Found complex',prev_test)
        if ((not prev_test or not prev_test_recent(prev_test,Prev,row)) and
            ('species' in row[group_field] or 'complex' in row[group_field])):
            prev_test = find_prev_for_generic(row['t1.patient_id'],row[group_field],Prev)
            if debug and row['t1.test_id'] == deb_test:
                print('Found species or generic',prev_test)
        if ((not prev_test or not prev_test_recent(prev_test,Prev,row)) and
            row['t4.class_broad'] in ['Mould','Dimorphic fungus','Other'] and
            'Unspecified organism' in Prev[row['t1.patient_id']]):
            prev_test = (row['t1.patient_id'],'Unspecified organism')
            if debug and row['t1.test_id'] == deb_test:
                print('Found unspecified organism',prev_test)
                print('Prev',Prev[row['t1.patient_id']])
                print('Is prev test recent?',prev_test_recent(prev_test,Prev,row))
        if ((not prev_test or not prev_test_recent(prev_test,Prev,row)) and
            row['t4.class_broad'] == 'Yeast' and
            row['t4.class_specific'] == 'Cryptococcus and related species' and
            'Unspecified organism' in Prev[row['t1.patient_id']]):
            prev_test = (row['t1.patient_id'],'Unspecified organism')
            #DEBUG
            #print('Found unspecified organism for cryptococcus',prev_test)
            if debug and row['t1.test_id'] == deb_test:
                print('Found unspecified organism',prev_test)
                print('Prev',Prev[row['t1.patient_id']])
                print('Is prev test recent?',prev_test_recent(prev_test,Prev,row))

        if ((not prev_test or not prev_test_recent(prev_test,Prev,row)) and
            row['t4.class_broad'] == 'Yeast' and
            row['t4.class_specific'] != 'Cryptococcus and related species' and
            'Unspecified organism' in Prev[row['t1.patient_id']]):
            # SREQ-120/SREQ-121: non-Cryptococcus Yeast (e.g. Meyerozyma) arriving after a
            # non-imaging "Unspecified organism" test should merge into that episode.
            # Time window is 90 days (Yeast), enforced by prev_test_recent() using the
            # incoming row's class_broad.
            # SREQ-114: do NOT link non-Cryptococcus Yeast to an 'Unspecified organism'
            # episode that was created by an imaging test (CT order_type=14 or
            # PET order_type=15).  Yeasts (excluding Cryptococcus) must never merge
            # with positive imaging tests.
            # Use has_imaging flag rather than order_type_id alone: a subsequent
            # non-imaging Unspecified organism test may have overwritten order_type_id
            # to a non-imaging value while the episode is still imaging-anchored.
            _unspec_entry = Prev[row['t1.patient_id']]['Unspecified organism']
            _unspec_has_imaging = _unspec_entry.get('has_imaging', _unspec_entry.get('order_type_id', '') in ('14', '15'))
            if not _unspec_has_imaging:
                prev_test = (row['t1.patient_id'],'Unspecified organism')

        for uns in ['Yeast','Mould']:
            if (not prev_test or not prev_test_recent(prev_test,Prev,row)) and 'Unspecified '+uns.lower() in Prev[row['t1.patient_id']] and row['t4.class_broad'] == uns:
                prev_test = (row['t1.patient_id'],'Unspecified '+uns.lower())
                if debug and row['t1.test_id'] == deb_test:
                    print('Found unspecified yeast or mould',prev_test)

        if (not prev_test or not prev_test_recent(prev_test,Prev,row)) and row[group_field] == 'Nakaseomyces glabratus' and 'Candida glabrata complex' in Prev[row['t1.patient_id']]:
            prev_test = (row['t1.patient_id'],'Candida glabrata complex')
            if debug and row['t1.test_id'] == deb_test:
                print('Found candida glabrata complex',prev_test)

        if not prev_test or not prev_test_recent(prev_test,Prev,row):
            return(None)
    return(prev_test)

def create_episode_group_org(merged_f,epi_org_f,group_field='t3.organism_id',id_field='episode_organism_id',cutoff_date=None):

    data = read_csv_to_dicts(merged_f)
    sorted_data = sort_dicts_by_id_and_date(data)
    Prev = {}
    epi_id = 1
    out_l= []
    date_format='%Y-%m-%d'
    
    # Track current patient for debug output
    current_patient = None
    patient_tests = []
    
    for row in sorted_data:
        # Debug: Track patient change and output summary
        debug = False
        deb_pat = '391'
        
        if debug and current_patient and current_patient != row['t1.patient_id']:
            # Patient changed - output summary of previous patient
            if current_patient == deb_pat:
                print(f"\n=== Patient {current_patient} - Test Processing Summary ===")
                for test_info in patient_tests:
                    print(f"  Test {test_info['test_id']:>6} | Date: {test_info['date']} | Species: {test_info['species']:30} | {'LINKED to ' + test_info['linked_to'] if test_info['linked_to'] else 'NEW EPISODE':40} | Epi_org: {test_info['epi_org']}")
                print(f"=== End Patient {current_patient} ===\n")
            patient_tests = []
        
        current_patient = row['t1.patient_id']
        
        #DEBUG
        #if row['t1.patient_id'] == '5341':
            #print(row['t1.test_id'],row['t1.collection_date'])
            #print(row)
        # Skip Colonisation tests - not clinically significant for episode grouping
        if row['system_ifi_classification'] == 'Colonisation':
            continue
        date_row = datetime.strptime(row['t1.collection_date'],date_format)

        # Episodes dated before the retention cutoff are being purged from
        # the DB and must not be recreated on a later run. Excluding these
        # tests here (before any grouping) means every episode this pass can
        # produce is anchored no earlier than the cutoff, by construction.
        if cutoff_date and date_row < cutoff_date:
            continue

        # Debug for test_id 112562 - show what this test looks like
        #if row['t1.test_id'] == '112562':
        #    print(f"DEBUG 112562: Processing test - species={row['t4.species']}, class_broad={row['t4.class_broad']}, order_type_id={row.get('t1.order_type_id', 'missing')}, collection_date={row['t1.collection_date']}")
        
        prev_test = get_prev_test(row,group_field,Prev)

        if prev_test:             
             (epi_org,epi_id,row,Prev) = get_epi_org_and_prev(epi_id,epi_org,prev_test,Prev,row,date_row)
             # Debug: Record linked test
             if debug and current_patient == deb_pat:
                 patient_tests.append({
                     'test_id': row['t1.test_id'],
                     'date': row['t1.collection_date'],
                     'species': row['t4.species'],
                     'linked_to': f"test with species '{prev_test[1]}'",
                     'epi_org': epi_org
                 })
        else:
            epi_org = epi_id
            epi_id += 1
            Daux = {'date':date_row,'first_date':date_row,'epi_org':epi_org,'class_broad':row['t4.class_broad'],'class_specific':row['t4.class_specific'],'final_infection_site_id':row['final_infection_site_id'],'order_type_id':row['t1.order_type_id']}
            if row['t1.patient_id'] in Prev:
                Prev[row['t1.patient_id']][row[group_field]] = Daux
            else:
                Prev[row['t1.patient_id']] = {row[group_field]:Daux}
            # Debug: Record new episode
            if debug and current_patient == deb_pat:
                patient_tests.append({
                    'test_id': row['t1.test_id'],
                    'date': row['t1.collection_date'],
                    'species': row['t4.species'],
                    'linked_to': None,
                    'epi_org': epi_org
                })
        row[id_field] = epi_org
        row['system_classification_id'] = row['system_ifi_classification_id']
        row['system_diagnosis_date'] = row['t1.collection_date']
        row['ifi_organism_site_id'] = row['final_infection_site_id']
        
        del row['t1.collection_date']
        del row['system_ifi_classification_id']
        del row['final_infection_site_id']
        del row['t1.order_id']

        del_rows = ['t1.encounter_id', 't1.is_sterile', 't6.specimen_source', 't1.order_type_id', 't7.order_name', 't3.organism_id', 't5.specimen_type', 't1.order_name_id']
        for r in del_rows:
            del row[r]
        out_l.append(row)
    
    # Debug: Output final patient summary if we ended on the debug patient
    debug = False
    deb_pat = '391'
    if debug and current_patient == deb_pat and patient_tests:
        print(f"\n=== Patient {current_patient} - Test Processing Summary ===")
        for test_info in patient_tests:
            print(f"  Test {test_info['test_id']:>6} | Date: {test_info['date']} | Species: {test_info['species']:30} | {'LINKED to ' + test_info['linked_to'] if test_info['linked_to'] else 'NEW EPISODE':40} | Epi_org: {test_info['epi_org']}")
        print(f"=== End Patient {current_patient} ===\n")
    
    out_l = sorted(out_l, key=lambda x: (int(x['t1.patient_id']), int(x['episode_organism_id'])))
    with open(epi_org_f, 'w', newline='', encoding='utf-8') as outfile:
        #h = ['episode_organism_id','t1.test_id','t1.encounter_id','t1.patient_id','t1.order_id','t3.organism_id','t4.species','t4.genus','t4.class_broad','t4.class_specific','t1.is_sterile','t6.specimen_source','t5.specimen_type','t1.order_name','t1.collection_date','system_ifi_classification_id','final_infection_site_id','infection_site_name','system_classification_id','system_classification_id_text']
        #h = [id_field,'system_diagnosis_date','t1.test_id','t1.encounter_id','t1.patient_id','t1.order_type_id','t3.organism_id','t4.species','t4.genus','t4.class_broad','t4.class_specific','t1.is_sterile','t6.specimen_source','t5.specimen_type','t7.order_name','t1.order_name_id','ifi_organism_site_id','infection_site_name','system_classification_id','system_ifi_classification']
        h = ['t1.patient_id','t1.test_id','system_diagnosis_date','infection_site_name','t4.species','t4.genus','t4.class_specific','t4.class_broad','system_ifi_classification','ifi_organism_site_id',id_field,'system_classification_id']
        writer = csv.DictWriter(outfile, fieldnames=h)
        writer.writeheader()
        for r in out_l:
            #if 'system_diagnosis_date' in r:
            #    r['system_diagnosis_date'] = format_date_dmy(r['system_diagnosis_date'])
            writer.writerow(r)

def update_site(epi_out_site,e,group_id):
    site = e['ifi_organism_site_id']
    site_n = e['infection_site_name']
    id = e[group_id]
    if id not in epi_out_site.keys():
        epi_out_site[id] = {}
        epi_out_site[id]['ifi_organism_site_id'] = [site]
        epi_out_site[id]['infection_site_name'] = [site_n]
    else:
        if site not in epi_out_site[id]['ifi_organism_site_id']:
            epi_out_site[id]['ifi_organism_site_id'].append(site)
            epi_out_site[id]['infection_site_name'].append(site_n)
    return(epi_out_site)

def update_date(epi_date,e,group_id):
    edate = e['system_diagnosis_date']
    if e[group_id] in epi_date:
        if edate < epi_date[e[group_id]]:
            epi_date[e[group_id]] = edate
    else:
        epi_date[e[group_id]] = edate
    return(epi_date)

def update_system_classification_id(epi_class, e, group_id, test_order_type_ids=None):
    # SREQ-127: Imaging tests (CT order_type=14, PET order_type=15) contribute Possible
    # classification but should not override the classification earned by microbiological
    # tests (e.g. Proven from a culture/PCR).  Skip imaging tests here; they are still
    # collected in epi_test and used by check_possible_to_probable_upgrade.
    if test_order_type_ids is not None:
        order_type = test_order_type_ids.get(e['t1.test_id'])
        if order_type in ('14', '15'):
            return epi_class

    class_id = (e['system_classification_id'], e['system_ifi_classification'])
    if e[group_id] in epi_class:
        if class_id[0] < epi_class[e[group_id]][0]:
            epi_class[e[group_id]] = class_id
    else:
        epi_class[e[group_id]] = class_id
    return epi_class

def check_possible_to_probable_upgrade(test_ids, test_order_type_ids):
    """
    Check if episode with 'Possible' classification should be upgraded to 'Probable'.
    Requires: at least one CT test (order_type == 14) + at least one of [2,7,8,11]
    """
    has_ct = False
    has_required_test = False
    required_test_types = {'2', '7', '8', '9','11'}
    
    for test_id in test_ids:
        if test_id in test_order_type_ids:
            order_type = test_order_type_ids[test_id]
            if order_type == '14':
                has_ct = True
            elif order_type in required_test_types:
                has_required_test = True
    
    return has_ct and has_required_test

def create_episode_group_singlerow(in_csv,out_csv,group_id,test_order_type_ids=None):
    epi = read_csv_to_dicts(in_csv)
    
    epi_out = {}
    epi_out_site = {}
    epi_test = {}
    epi_date = {}
    epi_class = {}
    for e in epi:
        epi_class = update_system_classification_id(epi_class, e, group_id, test_order_type_ids)
        epi_date = update_date(epi_date,e,group_id)
        epi_out_site = update_site(epi_out_site,e,group_id)
        if e[group_id] in epi_test:
            epi_test[e[group_id]].append(e['t1.test_id'])
        else:
            epi_test[e[group_id]] = [e['t1.test_id']]

        if e[group_id] in epi_out:
            # Compare against the species already stored for THIS group, not a
            # bare loop variable -- a bare variable only reflects the previous
            # row of the same group when the input happens to be pre-sorted by
            # group_id (true for this function's one current caller, but not
            # enforced or checked here, so it silently broke for any other
            # caller/order). Reading it back from epi_out is correct regardless
            # of input order.
            current_species = epi_out[e[group_id]]['row']['t4.species']
            if e['t4.species'] == 'Unspecified organism':
                continue
            if 'Unspecified ' in e['t4.species']:
                if 'Unspecified organism' not in current_species:
                    continue
            if ' species' in e['t4.species']:
                if 'Unspecified ' not in current_species:
                    continue
            if ' complex' in e['t4.species']:
                if 'Unspecified ' not in current_species and ' species' not in current_species:
                    continue

        epi_out[e[group_id]] = {'row':e}

    with open(out_csv, 'w', newline='', encoding='utf-8') as outfile:
        #h = [group_id,'system_diagnosis_date','ifi_organism_site_id','infection_site_name','t1.test_id','t1.encounter_id','t1.patient_id','t3.organism_id','t4.species','t4.genus','t4.class_broad','t4.class_specific','t1.is_sterile','t6.specimen_source','t5.specimen_type','t1.order_name','t7.order_name', 't1.order_name_id', 't1.order_type_id','system_classification_id','system_ifi_classification']
        h = ['t1.patient_id',group_id,'t1.test_id','system_diagnosis_date','infection_site_name','t4.species','t4.genus','t4.class_specific','t4.class_broad','system_classification_id','system_ifi_classification']
        #h = [group_id,'system_diagnosis_date','t1.test_id','t1.patient_id','t4.class_broad','t6.specimen_source','t5.specimen_type','t1.order_name','ifi_organism_site_id','infection_site_name','system_classification_id','system_classification_id_text']
        writer = csv.DictWriter(outfile, fieldnames=h)
        writer.writeheader()
        for e in epi_out:
            row = epi_out[e]['row']
            row['ifi_organism_site_id'] = ';'.join(epi_out_site[e]['ifi_organism_site_id'])
            # Remove 'Unknown' if there are multiple infection_site_name terms
            site_names = epi_out_site[e]['infection_site_name']
            if len(site_names) > 1:
                site_names = [s for s in site_names if s.strip().lower() != 'unknown']
            row['infection_site_name'] = ';'.join(sorted(site_names))
            row['t1.test_id'] = ';'.join(epi_test[e])
            row['system_diagnosis_date'] = epi_date[e]
            # Imaging-only episode_organism_ids have no entry in epi_class (all tests were
            # skipped as order_type 14/15). Fall back to Possible so SREQ-128 can still
            # upgrade it if the required non-imaging test is present in the same group.
            if e in epi_class:
                row['system_classification_id'] = epi_class[e][0]
                row['system_ifi_classification'] = epi_class[e][1]
            else:
                row['system_classification_id'] = '1'
                row['system_ifi_classification'] = 'Possible'
            
            # Check if 'Possible' or 'Pending' should be upgraded to 'Probable'.
            # Requires CT (order_type 14) + at least one non-imaging supporting
            # test (order_type in {2,7,8,9,11}).  Pending is included because the
            # XTC-310 imaging-skip causes groups whose only Possible evidence was
            # a CT scan to fall back to Pending (the classification of the remaining
            # non-imaging tests).  Proven and Colonisation are NOT overridden.
            if (test_order_type_ids and 
                row['system_ifi_classification'] in ('Possible', 'Pending') and
                check_possible_to_probable_upgrade(epi_test[e], test_order_type_ids)):
                #print(f"Upgrading patient {row['t1.patient_id']} episode {e} from '{row['system_ifi_classification']}' to 'Probable'")
                row['system_classification_id'] = '2'
                row['system_ifi_classification'] = 'Probable'
            
            #if 'system_diagnosis_date' in row:
            #    row['system_diagnosis_date'] = format_date_dmy(row['system_diagnosis_date'])
            del row['ifi_organism_site_id']
            writer.writerow(row)

def group_episode_org(outdir='Data/tmp', cutoff_date=None):
    """
    Group tests into episodes by organism, from the full, unfiltered test
    data. Review status (is_ifi_reviewer_agreed / clinical_review_status) is
    not considered here -- it is preserved per episode later, when
    update_episodes.py matches a grouping back to its existing DB episode.

    cutoff_date, if given, excludes tests collected before that date from
    grouping entirely, so no episode produced by this run can be anchored
    earlier than it -- used to stop episodes purged from the DB (retention
    cutoff) from being silently recreated on the next run.
    """
    input_file = f'{outdir}/tmp_merged2.csv'
    test_order_type_ids = get_test_order_type_ids(input_file)
    create_episode_group_org(input_file,f'{outdir}/epi_org.csv','t4.species','episode_organism_id',cutoff_date=cutoff_date)
    create_episode_group_singlerow(f'{outdir}/epi_org.csv',f'{outdir}/epi_org_table.csv','episode_organism_id',test_order_type_ids)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combine tests into episodes")
    # --password/--database/--server are no longer used by this script (episode
    # building no longer queries the DB), but are kept so full_pipeline_scheduler.py's
    # existing invocation (which always passes them) does not need to change.
    parser.add_argument("--password", type=str, help="Unused (kept for CLI compatibility)")
    parser.add_argument("--database", type=str, default="IFI_DATABASE_NAME", help="Unused (kept for CLI compatibility)")
    parser.add_argument("--server", type=str, default="DB_SERVER_HOSTNAME", help="Unused (kept for CLI compatibility)")
    parser.add_argument("--outdir", type=str, default="Data/tmp", help="Output directory for intermediate CSV files")
    parser.add_argument("--cutoff-date", type=str, default=None,
                        help="YYYY-MM-DD. Tests collected before this date are excluded from "
                             "episode grouping, so episodes anchored earlier than it are never "
                             "(re)created. Omit to group on the full, unfiltered test history.")
    args = parser.parse_args()

    cutoff_date = datetime.strptime(args.cutoff_date, '%Y-%m-%d') if args.cutoff_date else None

    group_episode_org(args.outdir, cutoff_date=cutoff_date)
    group_episode(args.outdir)
