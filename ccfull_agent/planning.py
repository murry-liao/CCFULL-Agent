"""Bounded scan proposals, finite-grid inference and read-only failure diagnostics."""
import copy
import itertools
import math
import re
from pathlib import Path
from .model import require, number, digest, read_json, validate, ValidationError
from .engine import make_plan, prepare, save, sha, parse_columns, check_approval
from .experiment import load_success, evaluate, check_dataset

def parameter_parent(config, path):
    require(isinstance(path,str) and re.fullmatch(
        r'(potential\.(V0_MeV|r0_fm|a_fm)|radial\.(Rmax_fm|dr_fm)|'
        r'couplings\.(projectile|target)\.modes\.[01]\.(energy_MeV|beta_C|beta_N|E2_MeV|beta2|beta4))',path),
        'Unsupported scan parameter; only existing continuous model parameters may vary')
    parts=path.split('.')
    value=config
    try:
        for part in parts[:-1]:
            value=value[int(part)] if isinstance(value,list) else value[part]
        require(parts[-1] in value, 'Scan parameter absent from base config')
    except (IndexError,KeyError) as exc:
        raise ValidationError('Scan path is absent from the base configuration') from exc
    return value,parts[-1]

def scan_configs(base, axes):
    validate(base)
    require(isinstance(axes,dict) and 1<=len(axes)<=3, 'Specify 1..3 scan axes')
    size=1
    for path, values in axes.items():
        parameter_parent(base,path)
        require(isinstance(values,list) and 1<=len(values)<=30, 'Each axis needs 1..30 values')
        for value in values:
            number(value,path)
        require(len(set(values))==len(values), 'Duplicate scan values')
        size*=len(values)
    require(size<=30, 'Scan exceeds 30 configurations; split and discuss separately')
    configs=[]
    for i, values in enumerate(itertools.product(*axes.values())):
        c=copy.deepcopy(base)
        c['id']=base['id']+f'__scan_{i:02d}'
        c['assumptions'].append('Explicit finite grid; no fitted physical conclusion or run authorization implied.')
        for path,value in zip(axes,values):
            obj,key=parameter_parent(c,path)
            obj[key]=value
        c['assumptions'].append('Trial overrides: '+', '.join(f'{p}={v}' for p,v in zip(axes,values))+
            '. These explicit values supersede matching inherited parameter statements; literature values are reference only. '
            'Unvaried beta_N remains fixed even when beta_C varies; equality is not imposed implicitly.')
        validate(c)
        configs.append(c)
    return configs

def prepare_scan(base, axes, out, purpose):
    configs=scan_configs(base,axes)
    result=prepare(make_plan(configs,purpose),out)
    save(Path(out)/'scan.json',{'schema_version':1,'base_config':base,'axes':axes,
         'plan_sha256':digest(read_json(Path(out)/'plan.json')),'status':'proposed_not_authorized'})
    return result

def fixed_config(config, axes):
    config=copy.deepcopy(config)
    config.pop('id')
    config.pop('assumptions')
    for path in axes:
        obj,key=parameter_parent(config,path)
        obj[key]=None
    return config

def fit_grid(paths, dataset, scan, out, interpolation='none', delta_chi2=1.0):
    """Select among a declared finite grid. No solver, continuous CI or posterior."""
    require(set(scan)=={'schema_version','base_config','axes','plan_sha256','status'} and scan['schema_version']==1,
            'Use scan.json from prepare-scan')
    check_dataset(dataset)
    require(dataset['metadata']['error_model']!='none', 'Grid fitting requires declared uncertainties/covariance')
    number(delta_chi2,'delta chi2 threshold',0)
    axes=scan['axes']
    configs=scan_configs(scan['base_config'],axes)
    expected={tuple(parameter_parent(c,p)[0][parameter_parent(c,p)[1]] for p in axes) for c in configs}
    require(1<=len(paths)<=len(expected), 'Unexpected number of runs')
    results=[load_success(p) for p in paths]
    require(len({r['binary_sha256'] for r in results})==1, 'Mixed solver builds')
    fixed=digest(fixed_config(scan['base_config'],axes))
    evaluated=[]
    seen=set()
    for path,result in zip(paths,results):
        require(digest(fixed_config(result['config'],axes))==fixed, 'Run changes an undeclared factor')
        point=tuple(parameter_parent(result['config'],p)[0][parameter_parent(result['config'],p)[1]] for p in axes)
        require(point in expected and point not in seen, 'Unexpected or duplicate grid point')
        # Exclude unapproved ad-hoc solver artifacts; benchmark receipts also cannot masquerade as research approval.
        authorization=read_json(Path(path)/'authorization.json')
        require(authorization.get('scope')=='exact-plan', 'Grid fits require approved research runs')
        parent=Path(path).parent
        plan=read_json(parent/'plan.json')
        receipt=read_json(parent/'approval.json')
        check_approval(plan,receipt)
        require(authorization==receipt, 'Run receipt differs from batch receipt')
        require(authorization.get('plan_sha256')==digest(plan), 'Run authorization/plan mismatch')
        require(any(digest(c)==result['config_sha256'] for c in plan['configs']), 'Config absent from approved plan')
        seen.add(point)
        score=evaluate(result,dataset,interpolation)
        evaluated.append({'parameters':dict(zip(axes,point)),'chi2':score['chi2'],
                          'run':str(path),'result_sha256':sha(Path(path)/'result.json')})
    evaluated.sort(key=lambda r:r['chi2'])
    best=evaluated[0]
    for r in evaluated:
        r['delta_chi2']=r['chi2']-best['chi2']
    selected=[r for r in evaluated if r['delta_chi2']<=delta_chi2]
    profile=[]
    for p,values in axes.items():
        profile.append({'parameter':p,'points':[{'value':v,'chi2_min':min(
            (r['chi2'] for r in evaluated if r['parameters'][p]==v),default=None)} for v in sorted(values)]})
    boundary=[p for p,values in axes.items() if len(values)>1 and best['parameters'][p] in (min(values),max(values))]
    report={'status':'finite_grid_evaluated','kind':dataset['metadata']['kind'],
        'dataset_sha256':dataset['data_sha256'],'scan_sha256':digest(scan),'interpolation':interpolation,
        'expected_points':len(expected),'evaluated_points':len(evaluated),
        'missing_points':[dict(zip(axes,p)) for p in sorted(expected-seen)],
        'best_sampled_point':best,'best_on_axis_boundary':boundary,'ranking':evaluated,'profiles':profile,
        'delta_chi2_threshold':delta_chi2,'sampled_support':selected,
        'limitations':['Best means only among supplied accepted runs; missing/failed points may change the result.',
                      'Sampled delta-chi2 support is not a continuous confidence interval or posterior.',
                      'No interpolation in parameter space, nuisance fit, physical-model discrepancy or degrees-of-freedom inference.',
                      'Boundary optimum or sparse grid requires discussion of a further scan, never automatic execution.']}
    root=Path(out)
    root.mkdir(parents=True,exist_ok=False)
    save(root/'fit.json',report)
    return report

def convergence_report(paths, out, rtol=0.01):
    from .analysis import compare_numeric
    number(rtol, 'relative tolerance', 0, 1)
    require(2<=len(paths)<=30, 'Convergence comparison needs 2..30 accepted runs')
    results=[load_success(p) for p in paths]
    axes={'radial.Rmax_fm':[], 'radial.dr_fm':[]}
    base=digest(fixed_config(results[0]['config'],axes))
    require(all(digest(fixed_config(r['config'],axes))==base for r in results), 'Convergence changes physical parameters or energy grid')
    require(len({r['binary_sha256'] for r in results})==1, 'Mixed solver builds')
    report={'reference':str(paths[0]), 'checks':[],
        'limitations':['Agreement on these radial grids does not establish channel-space or physical-model convergence.']}
    for path,r in zip(paths[1:],results[1:]):
        comparison=compare_numeric(results[0],r,rtol=rtol,atol=0,fields=('sigma_mb',),barrier_atol=None)
        report['checks'].append({'run':str(path),'radial':r['config']['radial'],**comparison})
    report['passed']=all(c['passed'] for c in report['checks'])
    from .experiment import new_json
    return new_json(out,report)

def diagnose_run(path):
    root=Path(path)
    status=read_json(root/'status.json')
    changes=[]
    for name,expected in status.get('artifacts',{}).items():
        require(Path(name).name==name, 'Invalid artifact path')
        if not (root/name).is_file() or sha(root/name)!=expected:
            changes.append(name)
    report={'run':str(root),'original_status':status['status'],'original_error':status.get('error'),
            'changed_artifacts':changes,'classification':[],'probability':None,
            'actions':[], 'promotes_failed_result':False}
    if changes:
        report['classification'].append('artifact_integrity_failure')
        report['actions'].append('Restore or investigate changed artifacts before numerical interpretation.')
        return report
    if (root/'s-wave.dat').exists():
        try:
            rows=parse_columns(root/'s-wave.dat')
            require(rows,'Empty probability file')
            bad=[{'row':i+1,'E_cm_MeV':e,'P':p} for i,(e,p) in enumerate(rows) if not -1e-8<=p<=1+1e-6]
            report['probability']={'minimum':min(p for e,p in rows),'maximum':max(p for e,p in rows),
                'invalid_count':len(bad),'examples':bad[:10],
                'note':'Rows include partial-wave integrations; no l column is available; do not call these all s waves.'}
            if bad:
                report['classification'].append('probability_bounds_failure')
                report['actions'].append('Propose independent dr/Rmax convergence checks; preserve thresholds and source. New runs require exact-plan confirmation.')
        except (ValidationError,ValueError) as exc:
            report['classification'].append('malformed_probability_output')
            report['actions'].append(str(exc))
    error=str(status.get('error','')).lower()
    if 'timeout' in error or 'timed out' in error:
        report['classification'].append('timeout')
    if status['status']!='succeeded' and not report['classification']:
        report['classification'].append('other_execution_or_validation_failure')
        report['actions'].append('Inspect stdout.log, stderr.log and the original validation error; return code alone is insufficient.')
    if not report['classification']:
        report['classification'].append('no_recorded_failure')
    return report
