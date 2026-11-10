"""Build the public research data for the site from the private research information system (RIS).

Usage: python3 site/build_research.py /path/to/ResearchInformationSystem

Writes, next to this script:
  research.json     papers in the RIS bibliography: metadata, themes, our summary, OpenAlex abstract and citations
  candidates.json   core candidates from the OpenAlex sweep (sweep/tiers.tsv, tier 1) that nobody has screened yet
  parameters.json   the parameter register (params/params.tsv)
  topics.json       the themes with their descriptions and counts
Only these fields leave the RIS. Internal fields stay behind: bib `note`, note sections other than Summary,
the claims ledger, reports and research notes. Prints a short summary."""
import csv
import json
import pathlib
import re
import sys

ris = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ris / 'scripts'))
from ris_lib import bib_entries, split_frontmatter  # noqa: E402

out = pathlib.Path(__file__).resolve().parent


def clean(s):
    s = re.sub(r'\\[a-zA-Z]+\{([^}]*)\}', r'\1', s or '')
    s = s.replace('{', '').replace('}', '').replace('\\&', '&').replace('\\%', '%').replace('--', '–').replace('**', '')
    return re.sub(r'\s+', ' ', s).strip()


def authors(a):
    names = [clean(n) for n in re.split(r'\s+and\s+', a or '') if n.strip()]
    short = []
    for n in names:
        short.append(n.split(',')[0].strip() if ',' in n else n.split()[-1] if n.split() else n)
    if len(short) > 3:
        return short[0] + ' et al.'
    return ', '.join(short)


def section(body, head):
    m = re.search(r'^## ' + re.escape(head) + r'\s*\n(.*?)(?=^## |\Z)', body, re.S | re.M)
    return re.sub(r'\s+', ' ', m.group(1)).strip() if m else ''


# themes from the vocabulary
vocab = (ris / 'VOCABULARY.md').read_text()
theme_desc = dict(re.findall(r'^\| `([a-z0-9-]+)` \| (.+?) \|$', vocab.split('## Kinds')[0], re.M))

# notes
notes = {}
for p in (ris / 'papers/notes').glob('*.md'):
    fm, body = split_frontmatter(p.read_text())
    if fm:
        notes[p.stem] = {'themes': fm.get('themes') or [], 'status': fm.get('status', ''),
                         'summary': section(body, 'Summary')}

# OpenAlex enrichment
oa = {}
if (ris / 'papers/openalex.tsv').exists():
    with (ris / 'papers/openalex.tsv').open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r.get('openalex_id'):
                oa[r['citekey']] = r
abstracts = {}
if (ris / 'papers/openalex_abstracts.jsonl').exists():
    for line in (ris / 'papers/openalex_abstracts.jsonl').open():
        d = json.loads(line)
        abstracts[d['citekey']] = d.get('abstract', '')

# sweep themes for bib entries the sweep found (tier 0)
sweep_rows = []
if (ris / 'sweep/tiers.tsv').exists():
    with (ris / 'sweep/tiers.tsv').open() as fh:
        sweep_rows = list(csv.DictReader(fh, delimiter='\t'))
sweep_by_doi = {r['doi']: r for r in sweep_rows if r['tier'] == '0' and r['doi']}

papers = []
bib = bib_entries((ris / 'papers/references.bib').read_text())
for key, f in bib.items():
    if 'false positive' in f.get('note', '').lower():
        continue
    n = notes.get(key, {})
    o = oa.get(key, {})
    doi = f.get('doi', '').lower()
    themes = list(n.get('themes') or [])
    if not themes and doi in sweep_by_doi:
        themes = sweep_by_doi[doi]['themes'].split(',')
    venue = f.get('journal') or f.get('booktitle') or f.get('institution') or f.get('publisher') or f.get('howpublished') or ''
    bibtex_fields = [(k, v) for k, v in f.items() if not k.startswith('_') and k not in ('note', 'kind', 'bias')]
    bibtex = '@%s{%s,\n%s\n}' % (f['_type'], key, ',\n'.join(f'  {k} = {{{v}}}' for k, v in bibtex_fields))
    try:
        year = int(f.get('year', '')[:4])
    except ValueError:
        year = None
    papers.append({
        'id': key, 'title': clean(f.get('title')), 'authors': authors(f.get('author')) or clean(f.get('institution')),
        'year': year, 'venue': clean(venue), 'kind': f.get('kind', ''), 'bias': f.get('bias', ''),
        'doi': doi, 'url': f.get('url') or (('https://doi.org/' + doi) if doi else ''),
        'themes': sorted(set(t for t in themes if t in theme_desc)), 'summary': n.get('summary', ''),
        'read': n.get('status', ''), 'abstract': abstracts.get(key, ''),
        'cited_by': int(o['cited_by_count']) if o.get('cited_by_count', '').isdigit() else None,
        'openalex': o.get('openalex_id', ''), 'oa_pdf': o.get('oa_pdf_url', ''), 'bibtex': bibtex,
    })
papers.sort(key=lambda p: (-(p['year'] or 0), p['title']))

screened = {}
if (ris / 'sweep/screen/decisions_all.tsv').exists():
    with (ris / 'sweep/screen/decisions_all.tsv').open() as fh:
        screened = {r['id']: r['decision'] for r in csv.DictReader(fh, delimiter='\t')}
candidates = []
for r in sweep_rows:
    if r['tier'] != '1' or screened.get(r['openalex_id']) in ('include', 'exclude'):
        continue
    candidates.append({'id': r['openalex_id'], 'title': clean(r['title']), 'year': int(r['year']) if r['year'] else None,
                       'venue': r['venue'], 'type': r['type'], 'cited_by': int(r['cited_by_count'] or 0),
                       'doi': r['doi'], 'themes': [t for t in r['themes'].split(',') if t in theme_desc],
                       'queries': int(r['n_queries'])})
candidates.sort(key=lambda c: (-c['queries'], -c['cited_by']))

params = []
with (ris / 'params/params.tsv').open() as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        note = re.sub(r'\s*\[batch [a-z_]+\]', '', r.get('note', '')).strip()
        params.append({k: r[k] for k in ('id', 'domain', 'parameter', 'scope', 'value', 'low', 'high', 'unit', 'year',
                                         'method', 'source_url', 'locator', 'verification', 'bias', 'status')} | {'note': note})

topics = []
for t, d in theme_desc.items():
    topics.append({'id': t, 'description': d, 'papers': sum(t in p['themes'] for p in papers),
                   'candidates': sum(t in c['themes'] for c in candidates)})

years = {}
for p in papers:
    if p['year'] and p['year'] >= 2006:
        years.setdefault(p['year'], [0, 0])[0] += 1
for c in candidates:
    if c['year'] and c['year'] >= 2006:
        years.setdefault(c['year'], [0, 0])[1] += 1
stats = {'papers': len(papers), 'candidates': len(candidates), 'parameters': len(params), 'topics': len(topics),
         'verified': sum(r['status'] == 'ok' for r in params),
         'years': {str(y): years[y] for y in sorted(years)},
         'sweep_works': len(sweep_rows), 'reserve': sum(r['tier'] == '2' for r in sweep_rows),
         'screened': len(screened), 'screened_in': sum(v == 'include' for v in screened.values()),
         'screened_out': sum(v == 'exclude' for v in screened.values())}
(out / 'stats.json').write_text(json.dumps(stats, separators=(',', ':')))

# the list page loads a small index; abstracts, summaries and BibTeX sit in 64 detail shards
DETAIL = ('summary', 'abstract', 'bibtex', 'url', 'doi', 'openalex', 'oa_pdf')
SHARDS = 64


def shard(key):
    return sum(ord(c) for c in key) % SHARDS   # same rule in research.html


shards = [{} for _ in range(SHARDS)]
index = []
for p in papers:
    shards[shard(p['id'])][p['id']] = {k: p[k] for k in DETAIL if p.get(k)}
    index.append({k: v for k, v in p.items() if k not in DETAIL} | {'has_summary': bool(p['summary']), 'has_abstract': bool(p['abstract'])})
detail_dir = out / 'research'
detail_dir.mkdir(exist_ok=True)
for i, s in enumerate(shards):
    (detail_dir / f'd{i:02d}.json').write_text(json.dumps(s, ensure_ascii=False, separators=(',', ':')))
full_papers = papers
papers = index

for name, data in (('research.json', papers), ('candidates.json', candidates), ('parameters.json', params), ('topics.json', topics)):
    (out / name).write_text(json.dumps(data, ensure_ascii=False, separators=(',', ':')))
print('=== Summary ===')
print(f'papers {len(full_papers)} (with themes {sum(bool(p["themes"]) for p in full_papers)}, abstracts {sum(bool(p["abstract"]) for p in full_papers)}, '
      f'summaries {sum(bool(p["summary"]) for p in full_papers)}) | candidates {len(candidates)} | parameters {len(params)} | topics {len(topics)}'
      f' | detail shards {SHARDS}, largest {max((detail_dir / f"d{i:02d}.json").stat().st_size for i in range(SHARDS)):,d} bytes')
for name in ('research.json', 'candidates.json', 'parameters.json', 'topics.json'):
    print(f'{name:16s} {(out / name).stat().st_size:>10,d} bytes')
