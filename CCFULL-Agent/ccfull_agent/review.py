"""Explicit retrospective acceptance of small C0 overshoot; original runs stay immutable."""
from pathlib import Path
from .model import require, read_json, digest, validate, text
from .engine import (SOURCE_SHA, check_approval, collect_run, parse_columns, save, sha,
                     verify_vendor, now)

POLICY = 'c0-cross-section-small-overshoot-v1'
UPPER_TOLERANCE = 2e-4
ORIGINAL_ERROR = 'Penetrability outside [0,1]: investigate numerical accuracy'

def reviewed_result(path):
    """Recheck all ordinary conditions, except the explicitly bounded upper P tolerance."""
    root=Path(path).resolve()
    verify_vendor()
    status=read_json(root/'status.json')
    require(status.get('status')=='failed' and status.get('error')==ORIGINAL_ERROR,
            'Review only applies to the recorded probability-bound failure')
    require(status.get('returncode')==0 and status.get('source_sha256')==SOURCE_SHA,
            'Review requires normal execution of the pinned source')
    required={'config.json','build.json','authorization.json','OUTPUT','cross.dat','spin.dat','s-wave.dat','stdout.log','stderr.log'}
    require(required<=set(status.get('artifacts',{})), 'Incomplete artifact manifest')
    for name,expected in status['artifacts'].items():
        require(Path(name).name==name and (root/name).is_file() and sha(root/name)==expected,
                f'Original artifact changed or missing: {name}')
    config=read_json(root/'config.json')
    validate(config)
    require(status['config_sha256']==digest(config), 'Original config hash mismatch')
    require(all(config['couplings'][side]['type']=='inert' for side in ('projectile','target')),
            'This exception is limited to uncoupled C0')
    plan=read_json(root.parent/'plan.json')
    receipt=read_json(root.parent/'approval.json')
    check_approval(plan,receipt)
    require(receipt==read_json(root/'authorization.json') and
            any(digest(c)==digest(config) for c in plan['configs']), 'Run is not in its approved plan')
    build=read_json(root/'build.json')
    require(build['identity']['source_sha256']==SOURCE_SHA, 'Build source mismatch')
    result=collect_run(root,config,probability_upper_tolerance=UPPER_TOLERANCE)
    penetration=parse_columns(root/'s-wave.dat')
    pmax=max(p for _,p in penetration)
    count=sum(p>1+1e-6 for _,p in penetration)
    require(count>0, 'No upper-bound overshoot requiring this exception')
    warning=(f'USER-ACCEPTED C0 NUMERICAL WARNING: Pmax={pmax:.16g}; '
             f'{count} rows exceed 1+1e-6. Original run failed strict validation. '
             'Cross-section comparison accepted; probabilities remain uncorrected and numerical issue unresolved.')
    result['warnings'].append(warning)
    result.update(config=config,config_sha256=digest(config),source_sha256=SOURCE_SHA,
                  binary_sha256=build['binary_sha256'],status='accepted_with_warning')
    result['review']={'policy':POLICY,'original_run':str(root),'original_status':'failed',
                      'P_max':pmax,'rows_above_strict_bound':count,
                      'scope':'cross_section_comparison_only','upper_tolerance':UPPER_TOLERANCE}
    provenance={'original_status_sha256':sha(root/'status.json'),
                'original_plan_sha256':digest(plan),'original_receipt_sha256':digest(receipt)}
    return result,provenance

def review_c0(paths, statement, actor, out):
    text(statement,'actual user acceptance statement')
    text(actor,'actor')
    require(1<=len(paths)<=30 and len({str(Path(p).resolve()) for p in paths})==len(paths),
            'Need 1..30 distinct original runs')
    checked=[reviewed_result(p) for p in paths]
    root=Path(out)
    root.mkdir(parents=True,exist_ok=False)
    reviews=[]
    for i,(result,provenance) in enumerate(checked):
        directory=root/f'{i:02d}'
        directory.mkdir()
        record={'schema_version':1,'policy':POLICY,'status':'accepted_with_warning',
                'actor':actor,'statement':statement,'recorded_at':now(),
                'original_run':result['review']['original_run'],**provenance,'result':result,
                'note':'Retrospective use decision, not a new calculation, identity authentication or numerical correction.'}
        record['review_sha256']=digest(record)
        save(directory/'review.json',record)
        reviews.append(str(directory))
    summary={'status':'accepted_with_warning','review_paths':reviews,'original_runs_unchanged':True,
             'scope':'cross_section_comparison_only','new_research_runs':0}
    save(root/'summary.json',summary)
    return summary

def load_reviewed(path):
    record=read_json(Path(path)/'review.json')
    require(record.get('schema_version')==1 and record.get('policy')==POLICY
            and record.get('status')=='accepted_with_warning','Unsupported review record')
    require(record.get('review_sha256')==digest({k:v for k,v in record.items() if k!='review_sha256'}),
            'Review changed after acceptance')
    text(record.get('statement'),'acceptance statement')
    text(record.get('actor'),'actor')
    result,provenance=reviewed_result(record['original_run'])
    require(all(record[k]==v for k,v in provenance.items()) and digest(result)==digest(record['result']),
            'Original run or review result changed; prior acceptance is stale')
    result['review']['record_sha256']=record['review_sha256']
    return result
