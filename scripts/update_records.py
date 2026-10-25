"""Regenerate the record lists from the data files and the live Covering Repository tables.

usage: python scripts/update_records.py [--work DIR] [--vkt FILE --vktm FILE]

  --work DIR   the local working folder holding lotto_new/ (SUMMARY.tsv, recipes/) and the component files.
               Every derived covering there whose size now equals the live table value (i.e. it has been
               uploaded) and that is not yet in data/derived.tsv is added: its recipe goes to derived/recipes/
               and its component list (sizes and sha1 hashes, not the files) to data/components.tsv.
  --vkt/--vktm use saved table files instead of downloading download_archive_vkt.aspx / download_archive.aspx.

Writes RECORDS.md, records/search_plain.md, records/search_lottery.md, records/derived.md, records/components.md.
Search records stay listed after someone beats them; their status column shows the current table value.
"""
import sys, os, json, hashlib, argparse, datetime, urllib.request, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, DER = os.path.join(ROOT, 'data'), os.path.join(ROOT, 'derived')
SITE = 'https://www.coveringrepository.com/'
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/130 Safari/537.36'}

ORIGINS = [  # prefix of the source path in the working folder -> where the component comes from
    ('keep/', 'own search covering (v1/v2)'),
    ('mkeep/', 'own search covering (wlsm2)'),
    ('lotto4/leaves/', 'La Jolla Covering Repository archive (Zenodo 19735294)'),
    ('data/mrec/', 'Covering Repository (downloaded record)'),
    ('data/rec/', 'Covering Repository (downloaded record)'),
    ('data/rec2/', 'Covering Repository (downloaded record)'),
    ('data/mrec2/', 'Covering Repository (downloaded record)'),
    ('reshrink/out/', 'own covering (re-shrunk with wlsm2big)'),
    ('lotto_new/', 'own derived covering'),
    ('lotto4/built', 'own derived covering'),
    ('construct', 'own derived covering (shrunk)'),
]


def parse_table(text):
    t = {}
    for l in text.splitlines():
        x = l.split()
        if len(x) >= 5 and x[4] != '-':
            t[tuple(map(int, x[:4]))] = int(x[4].replace('.', ''))
    return t


def live_tables(a):
    def get(path, page):
        if path: return open(path).read()
        return urllib.request.urlopen(urllib.request.Request(SITE + page, headers=UA), timeout=120).read().decode()
    return parse_table(get(a.vkt, 'download_archive_vkt.aspx')), parse_table(get(a.vktm, 'download_archive.aspx'))


def read_tsv(name):
    p = os.path.join(DATA, name)
    if not os.path.exists(p): return []
    lines = open(p, encoding='utf8').read().splitlines()
    head = lines[0].split('\t')
    return [dict(zip(head, l.split('\t'))) for l in lines[1:] if l.strip()]


def write_tsv(name, head, rows):
    with open(os.path.join(DATA, name), 'w', encoding='utf8', newline='\n') as f:
        f.write('\t'.join(head) + '\n')
        for r in rows: f.write('\t'.join(str(r.get(h, '')) for h in head) + '\n')


def leaves(n, acc):
    if n['op'] == 'leaf': acc.append(n)
    if 'child' in n: leaves(n['child'], acc)
    if n['op'] == 'product':
        for tm in n['terms']: leaves(tm[6], acc); leaves(tm[7], acc)
    return acc


def add_uploaded(work, live):
    rows = read_tsv('derived.tsv'); have = {(r['v'], r['k'], r['t'], r['m']) for r in rows}
    comps = {r['file']: r for r in read_tsv('components.tsv')}
    by_hash = {r['sha1']: r['file'] for r in comps.values()}
    os.makedirs(os.path.join(DER, 'recipes'), exist_ok=True)
    today = datetime.date.today().isoformat(); added = 0
    for l in open(os.path.join(work, 'lotto_new', 'SUMMARY.tsv'), encoding='utf8').read().splitlines()[1:]:
        r = l.split('\t')
        if r[0] not in ('lottery', 'plain'): continue
        v, k, t, m = map(int, r[1:5]); prev, ours = int(r[5]), int(r[6])
        if (str(v), str(k), str(t), str(m)) in have or live.get((v, k, t, m)) != ours: continue
        name = '%d,%d,%d,%d' % (v, k, t, m)
        tree = json.load(open(os.path.join(work, 'lotto_new', 'recipes', name + '.json')))
        store_recipe(work, name, tree, comps, by_hash)
        rows.append(dict(kind=r[0], v=v, k=k, t=t, m=m, previous=prev, ours=ours, construction=r[10], uploaded=today))
        added += 1
    save_lists(rows, comps)
    return added


def store_recipe(work, name, tree, comps, by_hash):
    """Record every leaf cover of a recipe (size, origin, sha1; the file itself is not published, since coverings
    downloaded from the Covering Repository are for personal use only) and save the recipe."""
    for leaf in leaves(tree, []):
        src = leaf['src'].replace('\\', '/'); data = open(os.path.join(work, src), 'rb').read()
        h = hashlib.sha1(data).hexdigest()
        if h not in by_hash:
            nb = sum(1 for x in data.splitlines() if x.strip())
            base = 'C_%d_%d_%d_%d_%d' % (leaf['v'], leaf['k'], leaf['t'], leaf['m'], nb); fn = base + '.txt'; i = 2
            while fn in comps: fn = '%s_%d.txt' % (base, i); i += 1
            origin = next((o for p, o in ORIGINS if src.startswith(p)), 'unknown')
            comps[fn] = dict(file=fn, v=leaf['v'], k=leaf['k'], t=leaf['t'], m=leaf['m'], blocks=nb, origin=origin, sha1=h)
            by_hash[h] = fn
        leaf['src'] = 'components/' + by_hash[h]
    json.dump(tree, open(os.path.join(DER, 'recipes', name + '.json'), 'w'), separators=(',', ':'))


def save_lists(rows, comps):
    key = lambda r: (int(r['t']), int(r['m']), int(r['v']), int(r['k']))
    write_tsv('derived.tsv', ['kind', 'v', 'k', 't', 'm', 'previous', 'ours', 'construction', 'uploaded'], sorted(rows, key=key))
    write_tsv('components.tsv', ['file', 'v', 'k', 't', 'm', 'blocks', 'origin', 'sha1'], sorted(comps.values(), key=key))


def describe(n):
    if n['op'] == 'product':
        v1 = n['v1']; v2 = n['v'] - v1
        return 'prod[%d|%d]: ' % (v1, v2) + ' + '.join('C(%d,%d,%d,%d)xC(%d,%d,%d,%d)' % (v1, a, b, c, v2, d, e, f) for a, b, c, d, e, f, _, _ in n['terms'])
    if 'child' in n:
        ch = n['child']; return '%s of C(%d,%d,%d,%d)' % (n['op'], ch['v'], ch['k'], ch['t'], ch['m'])
    return n['op']


def add_built(work, live, spec):
    """spec = builddir:targets[:shrunkdir], relative to work. A build dir holds v,k,t,m.txt + .json recipe per cover,
    targets lists 'v k t m live_before ...'. Every cover whose built or shrunk size is now the live value is added
    (or updated when the live value is smaller than what is listed)."""
    parts = spec.split(':'); bdir, tfile = parts[0], parts[1]; sdir = parts[2] if len(parts) > 2 else None
    rows = read_tsv('derived.tsv'); idx = {(r['v'], r['k'], r['t'], r['m']): r for r in rows}
    comps = {r['file']: r for r in read_tsv('components.tsv')}; by_hash = {r['sha1']: r['file'] for r in comps.values()}
    os.makedirs(os.path.join(DER, 'recipes'), exist_ok=True)
    nb = lambda p: sum(1 for x in open(p, 'rb') if x.strip())
    today = datetime.date.today().isoformat(); added = updated = 0
    for l in open(os.path.join(work, tfile)):
        x = l.split()
        if len(x) < 5: continue
        v, k, t, m, before = map(int, x[:5]); name = '%d,%d,%d,%d' % (v, k, t, m)
        bp = os.path.join(work, bdir, name + '.txt')
        if not os.path.exists(bp): continue
        built = nb(bp); sp = os.path.join(work, sdir, name + '.txt') if sdir else None
        shrunk = nb(sp) if sp and os.path.exists(sp) else None
        now = live.get((v, k, t, m))
        if now not in (built, shrunk): continue
        sk = (str(v), str(k), str(t), str(m)); old = idx.get(sk)
        if old and int(old['ours']) <= now: continue
        tree = json.load(open(os.path.join(work, bdir, name + '.json')))
        store_recipe(work, name, tree, comps, by_hash)
        con = describe(tree) + ('' if now == built else ' (recipe gives %d), then shrunk with wlsm2big' % built)
        row = dict(kind='plain' if t == m else 'lottery', v=v, k=k, t=t, m=m, previous=before, ours=now, construction=con, uploaded=today)
        if old: old.update({kk: str(vv) for kk, vv in row.items() if kk != 'previous'}); updated += 1
        else: rows.append(row); idx[sk] = row; added += 1
    save_lists(rows, comps)
    return added, updated


OWN_DERIVED = {}  # (v,k,t,m) -> size of our derived covering, filled in render()


def status(ours, now, key=None):
    if now is None: return 'not in table'
    if now == ours: return 'current record'
    if key is not None and OWN_DERIVED.get(key) == now: return 'improved by our derived covering (now %d)' % now
    return 'beaten (now %d)' % now if now < ours else 'table shows %d' % now


def table(head, rows):
    out = ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)]
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


def render(live, today):
    os.makedirs(os.path.join(ROOT, 'records'), exist_ok=True)
    W = lambda p, s: open(os.path.join(ROOT, p), 'w', encoding='utf8', newline='\n').write(s + '\n')
    cnt = {}
    for r in read_tsv('derived.tsv'): OWN_DERIVED[tuple(int(r[x]) for x in ('v', 'k', 't', 'm'))] = int(r['ours'])

    sp = read_tsv('search_plain.tsv'); rows = []; held = 0
    for r in sp:
        v, k, t, prev, ours = (int(r[x]) for x in ('v', 'k', 't', 'previous', 'ours'))
        st = status(ours, live.get((v, k, t, t)), (v, k, t, t)); held += st == 'current record'
        rows.append((v, k, t, prev, ours, prev - ours, r['method'], r['v1_value'], st))
    cnt['plain'] = (held, len(rows))
    W('records/search_plain.md', '# Plain coverings C(v,k,t) found by search\n\n'
      'Found by the weighted local search in [v1](../v1/) and [v2](../v2/), starting from the previous record covering. '
      '"Previous" is the table value before our upload. Where v2 improved a record that v1 had set, the v1 value is shown. '
      'Entries stay listed after someone else beats them; "status" shows the table on %s.\n\n' % today +
      table(['v', 'k', 't', 'previous', 'ours', 'gain', 'method', 'v1 value', 'status'], rows))

    sl = read_tsv('search_lottery.tsv'); rows = []; held = 0
    for r in sl:
        v, k, t, m, prev, ours = (int(r[x]) for x in ('v', 'k', 't', 'm', 'previous', 'ours'))
        st = status(ours, live.get((v, k, t, m)), (v, k, t, m)); held += st == 'current record'
        rows.append((v, k, t, m, prev, ours, prev - ours, r['method'], st))
    cnt['lottery'] = (held, len(rows))
    W('records/search_lottery.md', '# Lottery coverings C(v,k,t,m) found by search\n\n'
      'Found by [wlsm2](../vktm/) starting from the previous record covering, and verified exactly (CP-SAT infeasibility '
      'proof or exhaustive `mverify`). "Status" shows the table on %s.\n\n' % today +
      table(['v', 'k', 't', 'm', 'previous', 'ours', 'gain', 'method', 'status'], rows))

    comps = {r['file']: r for r in read_tsv('components.tsv')}
    used = collections.defaultdict(list); rows = []; held = 0
    for r in read_tsv('derived.tsv'):
        v, k, t, m, prev, ours = (int(r[x]) for x in ('v', 'k', 't', 'm', 'previous', 'ours'))
        name = '%d,%d,%d,%d' % (v, k, t, m)
        tree = json.load(open(os.path.join(DER, 'recipes', name + '.json')))
        files = sorted({l['src'].split('/')[-1] for l in leaves(tree, [])}, key=lambda f: [int(x) for x in f[2:-4].split('_')])
        for f in files: used[f].append(name)
        cl = ', '.join('C(%s,%s,%s,%s)=%s' % (comps[f]['v'], comps[f]['k'], comps[f]['t'], comps[f]['m'], comps[f]['blocks'])
                       for f in files)
        st = status(ours, live.get((v, k, t, m))); held += st == 'current record'
        rows.append((v, k, t, m, prev, ours, prev - ours, r['construction'], cl,
                     '[recipe](../derived/recipes/%s.json)' % name, st))
    cnt['derived'] = (held, len(rows))
    W('records/derived.md', '# Derived coverings\n\n'
      'Built from published coverings by the constructions described in [derived/](../derived/README.md): '
      'products over a split of the points, complement duality and one-point extensions. Each entry names every '
      'component covering used; every one is on the Covering Repository (at this size or better). The component files '
      'are not included here: download them from the Repository, check them against the sha1 hashes in '
      '[data/components.tsv](../data/components.tsv), then rebuild with `python derived/build.py v,k,t,m`. '
      'Plain coverings are the case m = t. "Status" shows the table on %s.\n\n' % today +
      table(['v', 'k', 't', 'm', 'previous', 'ours', 'gain', 'construction', 'components', 'recipe', 'status'], rows))

    rows = []
    for f, c in sorted(comps.items(), key=lambda x: [int(y) for y in x[0][2:-4].split('_')]):
        now = live.get(tuple(int(c[x]) for x in ('v', 'k', 't', 'm')))
        site = 'equal' if now == int(c['blocks']) else ('site has %d' % now if now is not None else 'not in table')
        rows.append((f, c['blocks'], c['origin'], site, len(used[f]),
                     ', '.join(used[f][:12]) + (' ...' if len(used[f]) > 12 else '')))
    W('records/components.md', '# Component coverings\n\n'
      'Every covering used as a building block by a derived covering, where it comes from, and how it compares with the '
      'Covering Repository table on %s ("site has N" means the site now holds a smaller covering, which can only make '
      'the derived covering smaller). sha1 hashes are in [data/components.tsv](../data/components.tsv).\n\n' % today +
      table(['file', 'blocks', 'origin', 'vs. table', 'used by', 'derived coverings'], rows))

    total = sum(h for h, _ in cnt.values())
    W('RECORDS.md', '# Records\n\n'
      'Coverings found with this code and accepted by the [Covering Repository](https://www.coveringrepository.com). '
      'On %s, **%d entries hold our covering**.\n\n' % (today, total) +
      table(['list', 'how they were found', 'current records', 'listed'], [
          ('[Plain coverings by search](records/search_plain.md)', 'weighted local search, v1/v2', *cnt['plain']),
          ('[Lottery coverings by search](records/search_lottery.md)', 'wlsm2', *cnt['lottery']),
          ('[Derived coverings](records/derived.md)', 'constructions from published coverings', *cnt['derived'])]) +
      '\n\nThe search records are original coverings. The derived coverings are built from component coverings that are '
      'all published on the Covering Repository; [records/components.md](records/components.md) lists each component '
      'and where it comes from. Lists are regenerated with `python scripts/update_records.py`.')
    return cnt


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--work'); ap.add_argument('--built', action='append', default=[]); ap.add_argument('--vkt'); ap.add_argument('--vktm')
    a = ap.parse_args()
    vkt, vktm = live_tables(a)
    live = dict(vktm); live.update({(v, k, t, t): b for (v, k, t, m), b in vkt.items() if m == t})
    if a.work: print('added', add_uploaded(a.work, live), 'derived coverings')
    for spec in a.built: print(spec, 'added/updated', add_built(a.work, live, spec))
    print(render(live, datetime.date.today().isoformat()))
