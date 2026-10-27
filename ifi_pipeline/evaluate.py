import argparse
import os
import pyodbc
import csv
from datetime import datetime, timedelta
from scipy.stats import norm

def wilson_ci(k, n, confidence=0.95):
    """Wilson score confidence interval for a binomial proportion k/n.

    Preferred over the normal approximation for proportions near 0 or 1
    and for the modest sample sizes typical of these evaluation runs.
    """
    if n == 0:
        return (0.0, 0.0)
    z = norm.ppf(1 - (1 - confidence) / 2)
    phat = k / n
    denom = 1 + z ** 2 / n
    centre = phat + z ** 2 / (2 * n)
    adj = z * ((phat * (1 - phat) / n + z ** 2 / (4 * n ** 2)) ** 0.5)
    lower = (centre - adj) / denom
    upper = (centre + adj) / denom
    return (max(0.0, lower), min(1.0, upper))

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

def convert_dates(rows, date_fields):
    for row in rows:
        for field in date_fields:
            if field in row:
                row[field] = format_date_dmy(row[field])
    return rows

def read_dict(patid,inf):
    outd = {}
    with open(inf, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        # If the file has a dedicated IFI_DATABASE_NAME patient ID column, use it
        # as the key and expose it as 'patient_id' so downstream code is consistent.
        surv_col = 'patient_id_DB_IFI_DATABASE_NAME'
        use_surv_col = patid == 'patient_id' and surv_col in (reader.fieldnames or [])
        for row in reader:
            if use_surv_col:
                key = row[surv_col]
                row['patient_id'] = key          # normalise to 'patient_id'
            else:
                key = row[patid]
            if key in outd:
                outd[key].append(row)
            else:
                outd[key] = [row]
    return(outd)

def close_dates(pred_date,gold_date):
    #dfg = "%m/%d/%Y"
    dfg = "%d/%m/%Y"
    dfp = "%d-%m-%Y"
    
    # Convert the date strings to datetime objects
    date1 = datetime.strptime(pred_date, dfp).date()
    date2 = datetime.strptime(gold_date, dfg).date()

    threshold = timedelta(days=30)
    difference = abs(date1 - date2)
    
    return difference <= threshold

def gold_found_within_f(Fields,field_last,gold_row,pred):
    all_tests = []
    found_tests = ""
    same_site = False
    same_classification = False
    found = ""
    found_pred = {}
    for field in Fields:
        for r in pred:

            all_tests += r['t1.test_id'].split(';')
            if r['t4.'+field] == gold_row[field]:
                pred_date = r['system_diagnosis_date']
                gold_date = gold_row['system_diagnosis_date']
                if close_dates(pred_date,gold_date):
                    found = field
                    found_pred = r
                    found_tests = r['t1.test_id']
                    if gold_row['infection_site'] == r['infection_site_name']:
                        same_site = True
                    if gold_row['ifi_classification'] == r['system_ifi_classification']:
                        same_classification = True
                    break
        if found or field == field_last:
            break
    
    # If no exact match found, check for unspecified fallback (matching the logic in pred_found_within_f)
    if not found:
        for r in pred:
            if 'Unspecified' in r['t4.class_broad']:
                pred_date = r['system_diagnosis_date']
                gold_date = gold_row['system_diagnosis_date']
                
                if close_dates(pred_date, gold_date):
                    found = 'Unspecified'
                    found_pred = r
                    found_tests = r['t1.test_id']
                    
                    # Check site and classification matches for unspecified cases
                    if gold_row['infection_site'] == r['infection_site_name']:
                        same_site = True
                    
                    if gold_row['ifi_classification'] == r['system_ifi_classification']:
                        same_classification = True
                    break
    
    # If still no match and gold has unspecified species, match any prediction with close dates
    if not found and 'Unspecified' in gold_row['species']:
        for r in pred:
            pred_date = r['system_diagnosis_date']
            gold_date = gold_row['system_diagnosis_date']
            
            if close_dates(pred_date, gold_date):
                found = 'Unspecified_Gold'
                found_pred = r
                found_tests = r['t1.test_id']
                
                # Check site and classification matches
                if gold_row['infection_site'] == r['infection_site_name']:
                    same_site = True
                
                if gold_row['ifi_classification'] == r['system_ifi_classification']:
                    same_classification = True
                break
    
    return(found,same_site,same_classification,found_tests,all_tests,found_pred)

def pred_found_within_f(Fields,field_last,pred_row,gold):
    
    found = ""
    pred_date = pred_row['system_diagnosis_date']
    for field in Fields:
        for gr in gold:
            if pred_row['t4.'+field] == gr[field]:
                gold_date = gr['system_diagnosis_date']
                if close_dates(pred_date,gold_date):
                    found = field
                    break
        if found or field == field_last:
            break
    if not found and 'Unspecified' in pred_row['t4.class_broad']:
        for gr in gold:
            gold_date = gr['system_diagnosis_date']
            if close_dates(pred_date,gold_date):
                found = 'Unspecified'
                #print('Found unspecified:',pred_row)
                break
    
    # If still no match, check if any gold has unspecified species (symmetric to gold_found_within_f)
    if not found:
        for gr in gold:
            if 'Unspecified' in gr['species']:
                gold_date = gr['system_diagnosis_date']
                if close_dates(pred_date,gold_date):
                    found = 'Unspecified_Gold'
                    break
    
    return(found)

def elem_found(field,gold_row,pred):
    found = False
    for r in pred:

        if r['t4.'+field] == gold_row[field]:
            pred_date = r['system_diagnosis_date']
            gold_date = gold_row['system_diagnosis_date']
            if close_dates(pred_date,gold_date):
                found = True
                break
    return(found)

def add_test_types(row,test_type_map,orders):
    if 'tests' in row:
        tlist = row['tests'].split(';')
    else:
        tlist = row['t1.test_id'].split(';')
    found = []
    for t in tlist:
        if t and t in test_type_map:
            found.append(test_type_map[t])
    for order in orders:
        if order in found:
            row[order] = 1
        else:
            row[order] = 0
    return(row)

def get_sensitivity(gold,pred,Fields,field,test_type_map,orders):
    TP = 0
    FN = 0
    same_class = 0
    diff_class = 0
    same_site = 0
    diff_site = 0
    tp_rows = []
    fn_rows = []
    used_preds = set()  # Track which predictions were used as TP
    
    for pat in gold:
        predpat = []
        if pat in pred:
            predpat = pred[pat]
        else:
            predpat = []
        
        for gold_row in gold[pat]:
            (found,ssite,sclass,found_tests,all_pat_tests,found_pred) = gold_found_within_f(Fields,field,gold_row,predpat) 
            
            if found:
                TP += 1
                gold_row['level'] = found
                gold_row['tests'] = found_tests
                gold_row['pred_system_diagnosis_date'] = found_pred['system_diagnosis_date']
                gold_row['pred_infection_site'] = found_pred['infection_site_name']
                gold_row['pred_ifi_classification'] = found_pred['system_ifi_classification']
                gold_row['pred_species'] = found_pred['t4.species']
                gold_row['pred_genus'] = found_pred['t4.genus']
                gold_row['pred_class_specific'] = found_pred['t4.class_specific']
                gold_row['pred_class_broad'] = found_pred['t4.class_broad']
                gold_row['pred_episode_organism_id'] = found_pred['episode_organism_id']
                gold_row = add_test_types(gold_row,test_type_map,orders)
                tp_rows.append(gold_row)
                if ssite:
                    same_site += 1
                else:
                    diff_site += 1
                if sclass:
                    same_class += 1
                else:
                    diff_class += 1
                predpat = [p for p in predpat if p != found_pred]
                used_preds.add(id(found_pred))  # Track this prediction as used
                #print('True positive:',gold_row)
            else:
                FN += 1
                gold_row['tests'] = ';'.join(all_pat_tests)
                gold_row = add_test_types(gold_row,test_type_map,orders)
                fn_rows.append(gold_row)
                #print('False negative:',gold_row)
                #print('Patient Pred:',predpat)

    return(TP,FN,same_class,diff_class,same_site,diff_site,tp_rows,fn_rows,used_preds)

def get_ppv(gold,pred,Fields,field,test_type_map,orders,used_preds):
    FP = 0
    fp_rows = []
    for pat in pred:
        for pred_row in pred[pat]:
            # Check if this prediction was used as a TP
            if id(pred_row) not in used_preds:
                FP += 1
                pred_row = add_test_types(pred_row,test_type_map,orders)
                fp_rows.append(pred_row)
    return(FP,fp_rows)

def date_filter(epis,dateformat):
    epis_out = {}
    for pat in epis:
        epi_out = []
        for epi in epis[pat]:
            date_object = datetime.strptime(epi['system_diagnosis_date'], dateformat)
            if 2018 <= date_object.year <= 2022:
                epi_out.append(epi)
        if epi_out:
            epis_out[pat] = epi_out
    return epis_out

def eval_sep_fields(gold,pred,Fields,test_type_map,orders,evaldir='Data/eval'):
    for field in Fields:
        (TP,FN,same_class,diff_class,same_site,diff_site,tp_rows,fn_rows,used_preds) = get_sensitivity(gold,pred,Fields,field,test_type_map,orders)
        (FP,fp_rows) = get_ppv(gold,pred,Fields,field,test_type_map,orders,used_preds)
        same_classp = same_class / (same_class + diff_class) if (same_class + diff_class) > 0 else 0
        same_sitep = same_site / (same_site + diff_site) if (same_site + diff_site) > 0 else 0
        print('\n'+field)
        print('TP:',TP)
        print(f'\tCorrect site %: {same_sitep:.2f}')
        print(f'\tCorrect ifi classification %: {same_classp:.2f}')
        print('FN:',FN)
        print('FP:',FP)
        rec = TP / (TP + FN)
        print(f'Sensitivity: {rec:.2f}')
        ppv = TP / (TP + FP)
        print(f'PPV: {ppv:.2f}')

        with open(f'{evaldir}/'+field+'_tp.csv', 'w', newline='') as csvfile:
            fieldnames = tp_rows[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(tp_rows)
        
        with open(f'{evaldir}/'+field+'_fn.csv', 'w', newline='') as csvfile:
            fieldnames = fn_rows[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(fn_rows)

        with open(f'{evaldir}/'+field+'_fp.csv', 'w', newline='') as csvfile:
            fieldnames = fp_rows[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(fp_rows)

def eval(gold,pred,Fields,test_type_map,orders,evaldir='Data/eval',confidence=0.95):
    os.makedirs(evaldir, exist_ok=True)
    (TP,FN,same_class,diff_class,same_site,diff_site,tp_rows,fn_rows,used_preds) = get_sensitivity(gold,pred,Fields,'class_broad',test_type_map,orders)
    (FP,fp_rows) = get_ppv(gold,pred,Fields,'class_broad',test_type_map,orders,used_preds)
    same_classp = same_class / (same_class + diff_class) if (same_class + diff_class) > 0 else 0
    same_sitep = same_site / (same_site + diff_site) if (same_site + diff_site) > 0 else 0
    same_class_lo, same_class_hi = wilson_ci(same_class, same_class + diff_class, confidence)
    same_site_lo, same_site_hi = wilson_ci(same_site, same_site + diff_site, confidence)

    with open(f'{evaldir}/result.txt', 'w') as result:

        ci_pct = int(confidence * 100)
        print('TP:',TP, file=result)
        print(f'\tCorrect site %: {same_sitep:.2f} ({ci_pct}% CI: {same_site_lo:.2f}-{same_site_hi:.2f})', file=result)
        print(f'\tCorrect ifi classification %: {same_classp:.2f} ({ci_pct}% CI: {same_class_lo:.2f}-{same_class_hi:.2f})', file=result)
        print('FN:',FN, file=result)
        print('FP:',FP, file=result)
        rec = TP / (TP + FN)
        rec_lo, rec_hi = wilson_ci(TP, TP + FN, confidence)
        print(f'Sensitivity: {rec:.2f} ({ci_pct}% CI: {rec_lo:.2f}-{rec_hi:.2f})', file=result)
        ppv = TP / (TP + FP)
        ppv_lo, ppv_hi = wilson_ci(TP, TP + FP, confidence)
        print(f'PPV: {ppv:.2f} ({ci_pct}% CI: {ppv_lo:.2f}-{ppv_hi:.2f})', file=result)

    date_fields = [
        'system_diagnosis_date',
        'pred_system_diagnosis_date'
    ]

    tp_rows = convert_dates(tp_rows, date_fields)
    fn_rows = convert_dates(fn_rows, date_fields)
    fp_rows = convert_dates(fp_rows, date_fields)

    with open(f'{evaldir}/tp.csv', 'w', newline='') as csvfile:
        fieldnames = ['gold_doc_id', 'patient_id', 'ifi_episode_id', 'if_episode_organism_id', 'system_diagnosis_date', 'pred_system_diagnosis_date','organism_id', 'species', 'pred_species', 'genus', 'pred_genus', 'class_specific', 'pred_class_specific', 'class_broad', 'pred_class_broad', 'ifi_infection_site_id', 'infection_site', 'pred_infection_site','system_ifi_classification_id', 'ifi_classification', 'pred_ifi_classification', 'pred_episode_organism_id', 'hospital', 'level', 'tests', 'Reason for exclusion', 'Antifungal Sensitivity', 'Aspergillus PCR', 'Biofire', 'CMV', 'CT scan', 'Cryptococcal antigen, Blood', 'Cryptococcal antigen, CSF', 'Culture', 'Full blood count', 'Galactomanan', 'Histopathology', 'Meningitis PCR', 'Other', 'PET scan', 'Panfungal PCR', 'TDM']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        tp_rows = sorted(tp_rows, key=lambda d: int(d['patient_id']))
        writer.writerows(tp_rows)

    with open(f'{evaldir}/fn.csv', 'w', newline='') as csvfile:
        fieldnames = fn_rows[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        fn_rows = sorted(fn_rows, key=lambda d: int(d['patient_id']))
        writer.writerows(fn_rows)

    with open(f'{evaldir}/fp.csv', 'w', newline='') as csvfile:
        fieldnames = fp_rows[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        fp_rows = sorted(fp_rows, key=lambda d: int(d['t1.patient_id']))
        writer.writerows(fp_rows)

def read_test_type_map(pwd, database='IFI_DATABASE_NAME', server='DB_SERVER_HOSTNAME'):
    # Database identity is site-specific and is read from the environment.
    # experiments/offline.py replaces this function with the exported CSVs, so
    # none of the experiments in this repository reach this code path.
    server_name = os.environ.get('IFI_DB_SERVER', server)
    db_name = os.environ.get('IFI_DB_DATABASE', database)
    uid = os.environ.get('IFI_DB_USER', 'DB_USERNAME')
    conn_str = (
        "DRIVER={ODBC Driver 17 for SQL Server};"
        f"SERVER={server_name};"
        f"DATABASE={db_name};"
        f"UID={uid};"
        f"PWD={pwd}"
    )

    conn = pyodbc.connect(conn_str)
    cursor = conn.cursor()
    cursor.execute("SELECT order_type,order_type_id FROM IFI.ref_Order_Type")
    rows = cursor.fetchall()
    maptype = {}
    for r in rows:
        maptype[str(r[1])] = r[0]

    cursor.execute("SELECT test_id,order_type_id FROM IFI.Test")
    rows = cursor.fetchall()
    conn.close()
    orders = []
    map = {}
    for r in rows:
        if r[1] is None:
            continue
        r[0] = str(r[0])
        r[1] = maptype[str(r[1])]
        map[r[0]] = r[1]
        if r[1] not in orders:
            orders.append(r[1])
    
    orders = sorted(orders)
    return(map,orders)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A simple example of using argparse")
    parser.add_argument("--password", type=str, nargs="?", help="The password")
    parser.add_argument("--database", type=str, default="IFI_DATABASE_NAME", help="Database name")
    parser.add_argument("--server", type=str, default="DB_SERVER_HOSTNAME", help="Database server name")
    parser.add_argument("--outdir", type=str, default="Data/tmp", help="Output directory for intermediate CSV files")
    parser.add_argument("--evaldir", type=str, default="Data/eval", help="Output directory for evaluation result files (default: Data/eval)")
    parser.add_argument("--confidence", type=float, default=0.95, help="Confidence level for Wilson score intervals on Sensitivity/PPV/site%%/classification%% (default: 0.95)")
    parser.add_argument("--goldfile", type=str,
                        default='gold_standard.csv',
                        help="Path to gold standard CSV (default: updateid file with IFI_DATABASE_NAME patient IDs)")
    args = parser.parse_args()

    if not args.password:
        print("Please provide a password using --password")

    gold = read_dict('patient_id', args.goldfile)
    # Strip automap-only columns and normalise organsim_id typo -> organism_id
    _automap_cols = {'_original_patient_id', '_shared_dates', '_overlap_pct_surv'}
    for pat in gold:
        for row in gold[pat]:
            for col in _automap_cols:
                row.pop(col, None)
            if 'organsim_id' in row and 'organism_id' not in row:
                row['organism_id'] = row.pop('organsim_id')
            # Strip trailing/leading whitespace from infection_site (gold standard data quality)
            if 'infection_site' in row:
                row['infection_site'] = row['infection_site'].strip()
    #gold = date_filter(gold,'%m/%d/%Y')
    gold = date_filter(gold,'%d/%m/%Y')
    pred = read_dict('t1.patient_id', f'{args.outdir}/epi_table.csv')
    #pred = date_filter(pred,'%Y-%m-%d')
    pred = date_filter(pred,'%d-%m-%Y')

    (test_type_map,orders) = read_test_type_map(args.password, args.database, args.server)
    Fields = ['species','genus','class_specific','class_broad']
    #eval_sep_fields(gold,pred,Fields,test_type_map,orders)
    eval(gold,pred,Fields,test_type_map,orders,evaldir=args.evaldir,confidence=args.confidence)

