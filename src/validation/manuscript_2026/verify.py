"""Independent Python checks of R models, survival curves and approved results."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm, chi2
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from common import REPOSITORY, require, write_json

EXPECTED = REPOSITORY/'data/manuscript_2026/expected_results.json'

def close(actual, expected, label, atol=1e-8, rtol=1e-7):
    require(np.allclose(actual, expected, atol=atol, rtol=rtol, equal_nan=True), f'Numerical mismatch: {label}')

def reml(beta, se):
    beta, variance = np.asarray(beta), np.square(se)
    def objective(tau2):
        weight = 1/(variance+tau2)
        mu = np.sum(weight*beta)/np.sum(weight)
        return .5*(np.log(variance+tau2).sum()+np.log(weight.sum())+np.sum(weight*(beta-mu)**2))
    maximum = max(4., float(np.var(beta)*10))
    optimum = minimize_scalar(objective, bounds=(0, maximum), method='bounded', options={'xatol':1e-14})
    require(optimum.success, 'Independent REML optimizer failed')
    tau2 = min([0., optimum.x], key=objective)
    weight = 1/(variance+tau2)
    beta_pool, se_pool = np.sum(weight*beta)/weight.sum(), np.sqrt(1/weight.sum())
    w0 = 1/variance
    typical_variance = (len(beta)-1)*w0.sum()/(w0.sum()**2-np.sum(w0**2))
    return dict(beta=beta_pool, se=se_pool, p=2*norm.sf(abs(beta_pool/se_pool)), tau2=tau2,
                I2=100*tau2/(tau2+typical_variance), weights=100*weight/weight.sum())

def survival_rows(frame):
    """KM and log-log Greenwood limits computed directly from event counts."""
    t, e = frame.time_months.to_numpy(), frame.event.to_numpy()
    probability, greenwood = 1., 0.
    rows=[]
    for time in np.unique(t):
        n = int(np.sum(t >= time)); events=int(np.sum((t==time)&(e==1))); censors=int(np.sum((t==time)&(e==0)))
        probability *= 1-events/n
        if events and n > events:
            greenwood += events/(n*(n-events))
        if 0 < probability < 1:
            loglog = np.log(-np.log(probability))
            delta = norm.ppf(.975)*np.sqrt(greenwood)/abs(np.log(probability))
            lower, upper = np.exp(-np.exp(loglog+delta)), np.exp(-np.exp(loglog-delta))
        elif probability == 1:
            lower = upper = 1.
        else:
            lower = upper = np.nan
        rows.append([time, probability, lower, upper, n, events, censors])
    return np.asarray(rows)

def logrank(frame):
    t, e, group = frame.time_months.to_numpy(), frame.event.to_numpy(), frame.high.to_numpy()
    difference, variance = 0., 0.
    for time in np.unique(t[e==1]):
        n = int(np.sum(t>=time)); n1 = int(np.sum((t>=time)&(group==1)))
        d = int(np.sum((t==time)&(e==1))); d1 = int(np.sum((t==time)&(e==1)&(group==1)))
        difference += d1 - d*n1/n
        if n > 1:
            variance += n1*(n-n1)*d*(n-d)/(n*n*(n-1))
    require(variance > 0, 'Non-estimable log-rank test')
    return difference*difference/variance

def fingerprint(frame):
    frame = frame.sort_values(['marker','group','time']).reset_index(drop=True)
    return hashlib.sha256(frame.to_csv(index=False, float_format='%.9f', na_rep='NA').encode()).hexdigest()

def verify(prepared, numerical, output):
    prepared, numerical, output = Path(prepared), Path(numerical), Path(output)
    expected = json.loads(EXPECTED.read_text())
    models = pd.read_csv(numerical/'model_statistics.csv').set_index('cohort')
    response = pd.read_csv(prepared/'response_analysis_data.csv')
    chemo = pd.read_csv(prepared/'chemo_analysis_data.csv')
    checks=[]
    for cohort, table in list(response.groupby('cohort'))+list(chemo.groupby('cohort')):
        q = np.nanquantile(table.fraction, [.25,.5,.75], method='linear')
        iqr = q[2]-q[0]
        require(iqr>0, 'Cannot verify per-IQR model with zero IQR')
        close(iqr, models.loc[cohort,'IQR'], cohort+' IQR')
        if cohort in set(response.cohort):
            x = (table.fraction-q[1])/iqr
            fit=sm.GLM(table.bad.to_numpy(), sm.add_constant(x.to_numpy()), family=sm.families.Binomial()).fit(maxiter=100, tol=1e-12)
            require(fit.converged, 'Independent logistic fit did not converge')
            beta, se = fit.params[1], fit.bse[1]
            unscaled=sm.GLM(table.bad.to_numpy(), sm.add_constant(table.fraction.to_numpy()), family=sm.families.Binomial()).fit(maxiter=100, tol=1e-12)
            close([beta,se], [unscaled.params[1]*iqr,unscaled.bse[1]*iqr], cohort+' raw scaling')
        else:
            table = table.dropna(subset=['OS_Time','OS_Status','fraction'])
            x = ((table.fraction-q[1])/iqr).to_numpy()[:,None]
            fit=PHReg(table.OS_Time.to_numpy(), x, status=table.OS_Status.to_numpy(), ties='efron').fit(disp=False)
            beta, se=fit.params[0],fit.bse[0]
        row=models.loc[cohort]
        close([beta,se], row[['beta','se']].astype(float), cohort+' independent coefficient/SE', atol=1e-7)
        close(np.exp(beta+norm.ppf([.025,.975])*se), row[['CI_low','CI_high']].astype(float), cohort+' interval', atol=1e-7)
        checks.append(dict(cohort=cohort, beta_difference=float(abs(beta-row.beta)), se_difference=float(abs(se-row.se))))

    for key in ['response_models','chemo_models']:
        for row in expected[key]:
            cohort=row['cohort'];actual=models.loc[cohort]
            for field in ['n','beta','se']:
                close(actual[field], row[field], cohort+' approved '+field)
            if key=='chemo_models':
                close(actual.joined_before_endpoint_omission, row['joined_before_endpoint_omission'], cohort+' pre-omission n')
                close(np.log(actual[['CI_low','CI_high']].astype(float)), [row['CI_low'],row['CI_high']], cohort+' historical CI')

    pools=pd.read_csv(numerical/'pooled_statistics.csv').set_index('pool')
    weights=pd.read_csv(numerical/'pool_weights.csv')
    pool_checks=[]
    for pool, roster, gold_key in [('chemo_OS',['monti','rchop','chop'],'chemo_pool'),('durability_primary',['moffitt','zuma7','stahl','li_cart'],'response_pool')]:
        m=models.loc[roster]; independent=reml(m.beta.to_numpy(),m.se.to_numpy());actual=pools.loc[pool]
        # Original chemotherapy metafor defaults stop earlier than the independent
        # bounded optimizer. Its historical output is checked more tightly below.
        tolerances = {'beta':2e-6,'se':2e-6,'p':2e-8,'tau2':1e-6,'I2':2e-3} if pool=='chemo_OS' else {field:1e-7 for field in ['beta','se','p','tau2','I2']}
        for field in ['beta','se','p','tau2','I2']:
            close(actual[field], independent[field], pool+' independent '+field, atol=tolerances[field], rtol=1e-7)
        for field in ['beta','se','p','tau2','I2','n']:
            close(actual[field],expected[gold_key][field],pool+' approved '+field, atol=1e-12 if field=='p' else 1e-8)
        close(actual.ratio, np.exp(actual.beta), pool+' ratio')
        close(actual[['CI_low','CI_high']].astype(float), np.exp(actual.beta+norm.ppf([.025,.975])*actual.se), pool+' Wald interval')
        gold_ci = [expected[gold_key]['CI_low'], expected[gold_key]['CI_high']]
        if pool=='chemo_OS':
            gold_ci = np.exp(gold_ci)  # Historical table stores log-HR limits.
        close(actual[['CI_low','CI_high']].astype(float), gold_ci, pool+' approved interval')
        w=weights[weights.pool.eq(pool)].set_index('cohort').loc[roster].weight_percent
        close(w,independent['weights'],pool+' independent weights',atol=3e-4 if pool=='chemo_OS' else 1e-7)
        exact_weights=1/(m.se**2+actual.tau2);exact_weights*=100/exact_weights.sum()
        close(w,exact_weights,pool+' weights')
        if pool=='durability_primary':
            close(w,[r['weight_percent'] for r in expected['response_models']],pool+' selected weights')
        pool_checks.append(dict(pool=pool, independent_beta=float(independent['beta']), independent_p=float(independent['p']), independent_I2=float(independent['I2']), historical_I2=float(actual.I2)))

    survival=pd.read_csv(prepared/'survival_analysis_data.csv')
    stats=pd.read_csv(numerical/'KM_statistics.csv').set_index('marker')
    curves=pd.read_csv(numerical/'km_curves.csv')
    risks=pd.read_csv(numerical/'risk_tables.csv')
    for marker, table in survival.groupby('marker'):
        row=stats.loc[marker]
        require(len(table)==22 and table.patient_id.nunique()==22,'LymphoMAP patient count')
        fit=PHReg(table.time_months.to_numpy(),table.high.to_numpy()[:,None],status=table.event.to_numpy(),ties='efron').fit(disp=False)
        close([fit.params[0],fit.bse[0]],row[['beta','se']].astype(float),marker+' independent KM Cox')
        statistic=logrank(table)
        close([statistic,chi2.sf(statistic,1)],row[['logrank_chisq','logrank_p']].astype(float),marker+' log-rank')
        for group, block in table.groupby('group'):
            c=curves[curves.marker.eq(marker)&curves.group.eq(group)]
            close(survival_rows(block),c[['time','survival','lower','upper','n_risk','n_event','n_censor']], marker+' '+group+' independent curve', atol=1e-8)
            r=risks[risks.marker.eq(marker)&risks.group.eq(group)]
            close([sum(block.time_months>=t) for t in r.time],r.n_risk,marker+' risk table')
    for row in expected['KM_statistics']:
        for field, value in row.items():
            if isinstance(value,(int,float)):
                close(stats.loc[row['marker'],field],value,row['marker']+' approved '+field)
    for name, table in [('km_curves',curves),('risk_tables',risks)]:
        require(len(table)==expected[name]['rows'] and fingerprint(table)==expected[name]['sha256'], 'Approved aggregate-table fingerprint mismatch: '+name)
    output.mkdir(parents=True, exist_ok=True)
    write_json(output/'independent_validation.json',dict(status='PASS',cohort_models=checks,pools=pool_checks,KM_curves_and_CIs=True,censor_counts=True,risk_tables=True,approved_KM_statistics=True,approved_curve_and_risk_fingerprints=True,reference_uncertainty_in_clinical_CI=False))
    print('Independent models, REML pools, KMs, censor counts and risk tables: PASS')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared-dir',type=Path,required=True)
    p.add_argument('--numerical-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();verify(args.prepared_dir,args.numerical_dir,args.output_dir)

if __name__=='__main__':
    main()
