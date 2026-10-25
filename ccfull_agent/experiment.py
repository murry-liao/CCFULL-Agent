"""Unit-explicit experimental imports and covariance-aware comparisons."""
import csv
import math
from pathlib import Path
import numpy as np
from .model import require, number, text, digest, read_json, validate, ValidationError
from .engine import save, sha, SOURCE_SHA

def new_json(path, data):
    require(not Path(path).exists(), 'Output already exists; choose a new destination')
    save(path, data)
    return data

def system_key(config):
    return {side: {k: config[side][k] for k in ('A', 'Z')} for side in ('projectile', 'target')}

def check_system(system):
    require(set(system)=={'projectile','target'}, 'System needs projectile and target')
    for side in system.values():
        require(set(side)=={'A','Z'}, 'Nucleus needs A and Z only')
        number(side['A'], 'A', 1, 300, True)
        number(side['Z'], 'Z', 1, side['A'], True)

def covariance_matrix(values, n):
    require(n<=1000, 'Full covariance currently limited to 1000 points')
    try:
        a = np.asarray(values, dtype=float)
        require(a.shape==(n,n) and np.isfinite(a).all(), 'Covariance must be finite N x N')
        require(np.allclose(a,a.T,rtol=1e-12,atol=0), 'Covariance must be symmetric')
        require(np.all(np.diag(a)>0), 'Covariance diagonal must be positive')
        np.linalg.cholesky(a)
    except (ValueError, np.linalg.LinAlgError) as exc:
        raise ValidationError('Covariance must be numeric and positive definite') from exc
    return a

def import_experiment(csv_path, metadata, out):
    required={'schema_version','id','kind','system','observable','source','energy_unit','energy_frame',
              'sigma_unit','uncertainty_convention','error_model'}
    require(set(metadata) in (required, required|{'covariance'}), 'Unexpected experimental metadata fields')
    require(metadata['schema_version']==1, 'Unsupported experimental schema')
    require(metadata['kind'] in ('measured','synthetic'), 'Explicit measured/synthetic label required')
    text(metadata['id'], 'dataset id')
    text(metadata['source'], 'data source with location')
    check_system(metadata['system'])
    require(metadata['observable']=='fusion', 'Fusion data required; capture/ER need an explicit physical conversion first')
    require(metadata['energy_unit'] in ('MeV','keV'), 'Energy unit must be MeV or keV')
    require(metadata['energy_frame'] in ('cm','lab_mass_number'), 'Use cm or explicit lab_mass_number approximation')
    require(metadata['sigma_unit'] in ('mb','b','microbarn'), 'Unsupported cross-section unit')
    require(metadata['error_model'] in ('none','independent','covariance'), 'Choose none, independent or covariance errors')
    require(metadata['uncertainty_convention']==('none' if metadata['error_model']=='none' else 'symmetric_1sigma'),
            'Only explicitly symmetric 1-sigma cross-section errors supported')
    factor={'mb':1.0,'b':1000.0,'microbarn':0.001}[metadata['sigma_unit']]
    ef=1 if metadata['energy_unit']=='MeV' else 0.001
    if metadata['energy_frame']=='lab_mass_number':
        ap, at = (metadata['system'][s]['A'] for s in ('projectile','target'))
        ef *= at/(ap+at)
    rows=[]
    with open(csv_path, newline='', encoding='utf-8-sig') as stream:
        reader=csv.DictReader(stream)
        expected=['energy','sigma'] + (['uncertainty'] if metadata['error_model']=='independent' else [])
        require(reader.fieldnames==expected, f'CSV columns must be exactly {expected}; energy errors/asymmetry are not silently ignored')
        for index, row in enumerate(reader):
            require(len(rows)<10000 and None not in row, 'Malformed CSV or >10000 points')
            try:
                e, s = float(row['energy'])*ef, float(row['sigma'])*factor
                u = float(row['uncertainty'])*factor if 'uncertainty' in row else None
            except (ValueError, TypeError) as exc:
                raise ValidationError(f'Invalid numeric CSV row {index+2}') from exc
            number(e, 'energy', 1e-12)
            number(s, 'cross section', 0)
            if u is not None:
                number(u, '1-sigma uncertainty', 1e-15)
            require(not rows or e>rows[-1]['E_cm_MeV'], 'Energies must be strictly increasing; duplicates require explicit treatment')
            rows.append({'E_cm_MeV':e, 'sigma_mb':s, 'uncertainty_mb':u})
    require(rows, 'Empty experimental dataset')
    cov=None
    if metadata['error_model']=='covariance':
        cov=(covariance_matrix(metadata.get('covariance'),len(rows))*factor**2).tolist()
        for i, row in enumerate(rows):
            row['uncertainty_mb']=math.sqrt(cov[i][i])
    else:
        require('covariance' not in metadata, 'Covariance supplied but not selected')
    data={'schema_version':1,'metadata':metadata,'rows':rows,'covariance_mb2':cov,
          'input_csv_sha256':sha(csv_path), 'conversion':{'energy_factor':ef,'sigma_factor':factor},
          'limitations':['No energy-uncertainty or asymmetric-error likelihood implemented.',
                          'No normalization or energy-shift nuisance parameters fitted.']}
    data['data_sha256']=digest(data)
    return new_json(out,data)

def check_dataset(dataset):
    require(dataset.get('data_sha256')==digest({k:v for k,v in dataset.items() if k!='data_sha256'}),
            'Dataset changed after import; re-import the source')
    check_system(dataset['metadata']['system'])
    require(dataset['rows'], 'Empty dataset')
    return dataset

def load_success(path):
    root=Path(path)
    status=read_json(root/'status.json')
    require(status['status']=='succeeded', f'Failed/incomplete run cannot be fitted: {path}')
    for name in ('result.json','config.json','build.json','OUTPUT','cross.dat','s-wave.dat','authorization.json'):
        require(name in status['artifacts'] and sha(root/name)==status['artifacts'][name], f'Run artifact changed: {name}')
    result=read_json(root/'result.json')
    validate(result['config'])
    require(result['status']=='succeeded' and result['source_sha256']==SOURCE_SHA, 'Unverified source or result status')
    require(digest(result['config'])==result['config_sha256']==status['config_sha256'], 'Config digest mismatch')
    require(result['config']==read_json(root/'config.json'), 'Result/config mismatch')
    build=read_json(root/'build.json')
    require(build['binary_sha256']==result['binary_sha256'] and build['identity']['source_sha256']==SOURCE_SHA,
            'Result/build identity mismatch')
    return result

def predictions(table, energies, interpolation='none'):
    require(interpolation in ('none','linear'), 'Interpolation must be none or linear')
    require(table, 'Empty solver table')
    x=np.array([r['E_cm_MeV'] for r in table],dtype=float)
    y=np.array([r['sigma_mb'] for r in table],dtype=float)
    require(np.isfinite(x).all() and np.isfinite(y).all() and np.all(y>=0) and np.all(np.diff(x)>0), 'Invalid solver grid')
    out=[]
    for e in energies:
        require(math.isfinite(e), 'Non-finite comparison energy')
        exact=np.where(np.isclose(x,e,rtol=0,atol=5.1e-6))[0]
        if len(exact)==1:
            out.append(float(y[exact[0]]))
        else:
            require(len(exact)==0 and interpolation=='linear', 'Different energy grids; explicitly choose linear interpolation')
            require(x[0]<=e<=x[-1], 'Extrapolation is forbidden')
            out.append(float(np.interp(e,x,y)))
    return out

def evaluate(result, dataset, interpolation='none'):
    check_dataset(dataset)
    require(system_key(result['config'])==dataset['metadata']['system'], 'Experimental and model systems differ')
    values=predictions(result['table'],[r['E_cm_MeV'] for r in dataset['rows']],interpolation)
    residual=np.array([p-r['sigma_mb'] for p,r in zip(values,dataset['rows'])])
    model=dataset['metadata']['error_model']
    chi2, whitened=None, None
    if model!='none':
        cov=(covariance_matrix(dataset['covariance_mb2'],len(values)) if model=='covariance'
             else np.diag([r['uncertainty_mb']**2 for r in dataset['rows']]))
        whitened=np.linalg.solve(np.linalg.cholesky(cov),residual)
        chi2=float(whitened@whitened)
        require(math.isfinite(chi2), 'Non-finite chi-squared')
    return {'id':result['config']['id'],'chi2':chi2,'n_points':len(values),'interpolation':interpolation,
            'residuals':[{'E_cm_MeV':r['E_cm_MeV'],'observed_mb':r['sigma_mb'],'predicted_mb':p,
                          'residual_mb':float(d),'uncertainty_mb':r['uncertainty_mb'],
                          'whitened_residual':None if whitened is None else float(whitened[i])}
                         for i,(r,p,d) in enumerate(zip(dataset['rows'],values,residual))],
            'interpretation':'Raw chi2 for supplied fixed predictions/covariance; no automatic degrees of freedom, p-value or mechanism selection.'}

def compare_experiment(paths, dataset, out, interpolation='none'):
    require(1<=len(paths)<=30, 'Compare 1..30 runs')
    results=[load_success(p) for p in paths]
    require(len({r['binary_sha256'] for r in results})==1, 'Mixed solver builds')
    scores=[evaluate(r,dataset,interpolation) for r in results]
    root=Path(out)
    root.mkdir(parents=True,exist_ok=False)
    report={'status':'analyzed','kind':dataset['metadata']['kind'],'dataset_sha256':dataset['data_sha256'],
            'runs':[{'path':str(p),'result_sha256':sha(Path(p)/'result.json')} for p in paths], 'comparisons':scores}
    save(root/'comparison.json',report)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes=plt.subplots(2,1,figsize=(7,7),sharex=True,layout='constrained')
    x=[r['E_cm_MeV'] for r in dataset['rows']]
    y=[r['sigma_mb'] for r in dataset['rows']]
    err=None if dataset['metadata']['error_model']=='none' else [r['uncertainty_mb'] for r in dataset['rows']]
    axes[0].errorbar(x,y,yerr=err,fmt='o',color='black',ms=4,label=dataset['metadata']['kind']+' data')
    for index,(result,score) in enumerate(zip(results,scores)):
        style=['-','--',':','-.'][index%4]
        axes[0].plot([r['E_cm_MeV'] for r in result['table']],[r['sigma_mb'] for r in result['table']],style,label=score['id'])
        axes[1].plot(x,[r['residual_mb'] for r in score['residuals']],style,label=score['id'])
        with (root/f'residuals-{index:02d}.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(score['residuals'][0]))
            writer.writeheader()
            writer.writerows(score['residuals'])
    axes[0].set_ylabel('Fusion cross section (mb)')
    axes[0].set_title(dataset['metadata']['kind'].upper()+' comparison — linear axes')
    axes[0].legend(fontsize=7)
    axes[1].axhline(0,color='grey',lw=.8)
    axes[1].set_ylabel('Model - data (mb)')
    axes[1].set_xlabel('Center-of-mass energy (MeV)')
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(root/'comparison.png',dpi=180)
    fig.savefig(root/'comparison.svg')
    plt.close(fig)
    save(root/'figure.json',{'kind':report['kind'],'dataset_sha256':dataset['data_sha256'],
         'errorbars':'Marginal symmetric 1-sigma; covariance correlations used in chi2, not visually represented.',
         'interpolation':interpolation,'scale':'linear; no points omitted','normalization':'none'})
    return {'report':str(root/'comparison.json'),'figure':str(root/'comparison.png'),'kind':report['kind']}
