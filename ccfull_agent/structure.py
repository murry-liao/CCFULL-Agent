"""Conservative NuDat ingestion and explainable screening, never solver execution."""
import math
import re
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import urlopen
from .model import require, number, beta_from_transition, digest
from .engine import save, sha, now
from .evidence import source_queries

SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
           'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce '
           'Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn '
           'Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()

def nucleus_info(name):
    m = re.fullmatch(r'(\d{1,3})([A-Z][a-z]?)', name)
    require(m is not None and m[2] in SYMBOLS, 'Use a real nucleus symbol, e.g. 144Sm')
    a, z = int(m[1]), SYMBOLS.index(m[2]) + 1
    require(z <= a <= 350, 'Invalid mass/charge combination')
    return {'name': name, 'A': a, 'Z': z}

def scalar(raw):
    """Only plain positive values and last-digit symmetric uncertainties are evaluated."""
    raw = ' '.join(raw.split())
    m = re.fullmatch(r'(\d+(?:\.\d*)?)(?:\s+(\d+))?', raw)
    result = {'raw': raw, 'value': None, 'uncertainty': None, 'status': 'unresolved'}
    if m:
        decimals = len(m[1].split('.')[1]) if '.' in m[1] else 0
        result.update(value=float(m[1]), uncertainty=int(m[2])*10**(-decimals) if m[2] else None,
                      status='parsed')
    return result

class Tables(HTMLParser):
    # NuDat mainTable omits some closing tr tags. Flush on opening tags too.
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables, self.table, self.row, self.cell = {}, None, [], None

    def flush_cell(self):
        if self.cell is not None:
            self.cell['text'] = ' '.join(self.cell['text'].split())
            self.row.append(self.cell)
            self.cell = None

    def flush_row(self):
        self.flush_cell()
        if self.table and self.row:
            self.tables[self.table].append(self.row)
        self.row = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'table':
            self.flush_row()
            self.table = attrs.get('id')
            if self.table:
                self.tables.setdefault(self.table, [])
        elif tag == 'tr':
            self.flush_row()
        elif tag in ('td', 'th') and self.table:
            self.flush_cell()
            self.cell = {'class': attrs.get('class', ''), 'text': '', 'rowspan': attrs.get('rowspan', '1')}
        elif self.cell is not None:
            if tag == 'br':
                self.cell['text'] += ' | '
            elif tag == 'i':
                self.cell['text'] += ' '
            elif tag == 'sup':
                # Keep scientific powers visibly unresolved, while preserving spin notation.
                if self.cell['text'].endswith('10'):
                    self.cell['text'] += '^'

    def handle_endtag(self, tag):
        if tag in ('td', 'th'):
            self.flush_cell()
        elif tag == 'tr':
            self.flush_row()
        elif tag == 'table':
            self.flush_row()
            self.table = None

    def handle_data(self, data):
        if self.cell is not None:
            self.cell['text'] += data

def cells(row, cls):
    return [c['text'] for c in row if cls in c['class'].split()]

def parse_nudat(path, nucleus, provenance):
    info = nucleus_info(nucleus)
    path = Path(path)
    require(path.stat().st_size <= 5_000_000, 'Snapshot exceeds 5 MB limit')
    require(provenance.get('sha256') == sha(path), 'Snapshot hash differs from provenance')
    require(provenance.get('url') and provenance.get('retrieved_at'), 'Source URL and retrieval date required')
    html = path.read_text(encoding='utf-8')
    require(re.search(r'List of levels for\s+' + re.escape(nucleus.upper()) + r'\s*</title>', html, re.I),
            'NuDat title does not identify the requested nucleus')
    parser = Tables()
    parser.feed(html)
    parser.flush_row()
    require('mainTable' in parser.tables and 'gammaTable' in parser.tables, 'Unsupported NuDat page layout')
    levels, branches, skipped = [], {}, []
    for index, row in enumerate(parser.tables['mainTable']):
        if any('header' in c['class'] for c in row):
            continue
        energies, spins = cells(row, 'elvl'), cells(row, 'jpi')
        if len(energies) != 1 or len(spins) != 1:
            skipped.append({'table': 'mainTable', 'row': index, 'reason': 'missing/ambiguous energy or spin'})
            continue
        energy = scalar(energies[0])
        level = {'energy_keV': energy, 'spin_parity': spins[0], 'transitions': [], 'source_row': index}
        levels.append(level)
        gammas, finals = cells(row, 'gamm'), cells(row, 'fin')
        if energy['value'] is not None and len(gammas) == 1 and len(finals) == 2:
            gs = [x.strip() for x in gammas[0].split('|') if x.strip()]
            fs = [x.strip() for x in finals[0].split('|') if x.strip()]
            if len(gs) == len(fs):
                branches[energy['value']] = [(scalar(g)['value'], scalar(f)['value']) for g, f in zip(gs, fs)]
    by_energy = {}
    for level in levels:
        by_energy.setdefault(level['energy_keV']['value'], []).append(level)
    initial, remaining = None, 0
    for index, row in enumerate(parser.tables['gammaTable']):
        if any('header' in c['class'] for c in row):
            continue
        starts = [c for c in row if 'gamm' in c['class'].split()]
        if starts:
            initial = scalar(starts[0]['text'])['value']
            require(starts[0]['rowspan'].isdigit(), 'Unsupported NuDat rowspan')
            remaining = int(starts[0]['rowspan'])
        if remaining <= 0:
            initial = None
        remaining -= 1
        extra = cells(row, 'extr')
        if not extra or 'B(E' not in extra[0]:
            continue
        gamma = [c['text'] for c in row if c['class'] == 'cell']
        candidates = by_energy.get(initial, []) if initial is not None else []
        if len(candidates) != 1 or len(gamma) != 1:
            skipped.append({'table': 'gammaTable', 'row': index, 'raw': extra[0], 'reason': 'ambiguous initial level/gamma'})
            continue
        g = scalar(gamma[0])
        # Match the tabulated gamma branch, not E(level)-E(gamma): recoil matters in light nuclei.
        destinations = [f for e, f in branches.get(initial, []) if e is not None and g['value'] is not None
                        and math.isclose(e, g['value'], rel_tol=0, abs_tol=0.00001)]
        for m in re.finditer(r'B\(E(\d+)\)\(W\.u\.\)([^,]*)', extra[0]):
            rhs = m[2].strip()
            strength = scalar(rhs[1:]) if rhs.startswith('=') else scalar(rhs)
            candidates[0]['transitions'].append({'lambda': int(m[1]), 'B_Wu': strength,
                'direction': 'down', 'gamma_keV': g, 'final_energy_keV': destinations[0] if len(destinations)==1 else None,
                'source_row': index, 'raw': m[0]})
    require(levels and any(l['energy_keV']['value'] == 0 for l in levels), 'No adopted ground state found')
    plain = unescape(re.sub('<[^>]+>', ' ', html))
    cutoff = re.search(r'Cutoff date:\s*([^|\n]+?)\s*(?:Full ENSDF|$)', plain)
    result = {'schema_version': 1, 'nucleus': info, 'source': dict(provenance),
            'evaluation_cutoff_raw': cutoff[1].strip() if cutoff else None,
            'levels': levels, 'unresolved_rows': skipped,
            'limitations': ['Only NuDat mainTable/gammaTable plain numeric fields are interpreted.',
                'Limits, asymmetric errors, tentative spins and complex expressions remain unresolved.',
                'No completeness claim for transition strengths; missing values are not zero.']}
    result['structure_sha256'] = digest(result)
    return result

def import_structure(path, nucleus, manifest, out):
    entries = [r for r in manifest if r['file'] == Path(path).name]
    require(len(entries) == 1, 'Need one matching snapshot manifest record')
    result = parse_nudat(path, nucleus, entries[0])
    require(not Path(out).exists(), 'Output already exists')
    save(out, result)
    return result

def fetch_structure(nucleus, out):
    nucleus_info(nucleus)
    root = Path(out)
    root.mkdir(parents=True, exist_ok=False)
    url = source_queries(nucleus)['nudat']
    save(root/'status.json', {'status': 'fetching', 'url': url})
    try:
        with urlopen(url, timeout=30) as response:
            data = response.read(5_000_001)
            require(len(data) <= 5_000_000, 'NuDat response exceeds 5 MB limit')
            require('text/html' in response.headers.get('Content-Type', ''), 'Expected HTML response')
            final_url = response.url
        path = root/'nudat.html'
        path.write_bytes(data)
        record = {'file': path.name, 'url': url, 'final_url': final_url, 'retrieved_at': now(), 'sha256': sha(path)}
        save(root/'manifest.json', [record])
        result = parse_nudat(path, nucleus, record)
        save(root/'structure.json', result)
        save(root/'status.json', {'status': 'succeeded'})
        return {'structure': str(root/'structure.json')}
    except Exception as exc:
        save(root/'status.json', {'status': 'failed', 'error': str(exc), 'url': url})
        raise

def beta_from_wu(B, A, Z, lam, radius, direction, uncertainty=None):
    """Single-particle W.u. (R_W=1.2 A^(1/3)); B direction supplied separately."""
    number(A, 'A', 1, 300, integer=True)
    number(lam, 'lambda', 2, 6, integer=True)
    unit = (3/(lam+3))**2 * (1.2*A**(1/3))**(2*lam) / (4*math.pi)
    number(B, 'B(W.u.)', 0)
    if uncertainty is not None:
        number(uncertainty, 'B uncertainty', 0)
    result = beta_from_transition(B*unit, 'e2_fm2lambda', A, Z, lam, radius, direction,
                                  None if uncertainty is None else uncertainty*unit)
    result.update(Wu_e2_fm2lambda=unit, weisskopf_radius_coefficient_fm=1.2,
                  convention='single-particle W.u.; direction converted once; coupling radius specified independently',
                  formula_source='https://www.nndc.bnl.gov/be2/help.html')
    return result

def assess_structure(structure, radius, max_energy=7.0, min_wu=5.0):
    require(structure.get('structure_sha256') == digest({k:v for k,v in structure.items() if k!='structure_sha256'}),
            'Structure changed after parsing; re-import the source snapshot')
    number(radius, 'coupling radius coefficient', 0.1, 5)
    number(max_energy, 'energy screen', 0.001, 100)
    number(min_wu, 'collectivity screening threshold', 0, 100000)
    info = structure['nucleus']
    require(nucleus_info(info['name']) == info, 'Nucleus metadata mismatch')
    levels = structure['levels']
    ground = [l for l in levels if l['energy_keV']['value'] == 0]
    require(len(ground) == 1, 'Ground-state assignment ambiguous')
    supported = ground[0]['spin_parity'] == '0+'
    selected = []
    for level in levels:
        e, spin = level['energy_keV']['value'], level['spin_parity']
        if e is None or e == 0 or e/1000 > max_energy or not any(s in spin for s in ('2+', '3-')):
            continue
        lam = {'2+': 2, '3-': 3}.get(spin)
        ts = [t for t in level['transitions'] if t['lambda']==lam and t['final_energy_keV']==0
              and t['B_Wu']['status']=='parsed' and t['B_Wu']['value']>0]
        gaps = ['Specify vibrational/rotational model with supporting literature.',
                'Specify beta_N, beta sign, phonon truncation and mutual-excitation choice.']
        beta, decision = None, 'needs_evidence'
        if not supported or lam is None:
            decision = 'blocked'
            gaps.insert(0, 'Pinned solver requires 0+ ground state and unambiguous multipolarity.')
        elif len(ts) == 1:
            b = ts[0]['B_Wu']
            beta = beta_from_wu(b['value'], info['A'], info['Z'], lam, radius, 'down', b['uncertainty'])
            decision = 'priority_candidate' if b['value'] >= min_wu else 'lower_priority_candidate'
        else:
            gaps.insert(0, 'Missing or conflicting resolved ground-state B(Eλ); no beta inferred.')
        selected.append({'energy_MeV': e/1000, 'spin_parity': spin, 'decision': decision,
                         'beta_C_magnitude': beta, 'evidence': ts, 'gaps': gaps})
    def first(spin):
        vals = [l['energy_keV']['value'] for l in levels if l['spin_parity']==spin and l['energy_keV']['value'] is not None and l['energy_keV']['value']>0]
        return min(vals) if vals else None
    e2, e4 = first('2+'), first('4+')
    return {'status': 'awaiting_physics_discussion', 'nucleus': info, 'source': structure['source'],
            'structure_sha256': structure['structure_sha256'],
            'solver_ground_state_supported': supported, 'candidates': selected,
            'R42': e4/e2 if e2 and e4 else None,
            'screening': {'max_energy_MeV': max_energy, 'min_downward_Wu': min_wu,
                          'meaning': 'User-adjustable triage heuristic, not proof a coupling is necessary; R42 alone never selects a rotor.'},
            'limitations': ['No beta_N, signed deformation, transfer strength or executable configuration is invented.',
                           'Magnitude uncertainty only propagates B at fixed coupling radius; model/radius uncertainty excluded.',
                           'Baseline inert and nested candidate comparisons require a separate reviewed plan.']}
