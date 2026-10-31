# Build the Deep Space Strategy site into docs/ for GitHub Pages: the landing page and the forwarding page that sends the
# old map addresses to solmaps (https://bhanneke.github.io/solmaps/), where the maps live since 2026-10-31.
import pathlib, shutil

root = pathlib.Path(__file__).resolve().parent.parent
docs = root / 'docs'
docs.mkdir(exist_ok=True)

def wrap(html):
    # pages are written without a document skeleton; standalone hosting needs one
    head = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"></head><body>\n')
    return head + html + '\n</body></html>\n'

(docs / 'index.html').write_text(wrap((root / 'site' / 'index.html').read_text()))
shutil.copy(root / 'site' / '404.html', docs / '404.html')
shutil.copy(root / 'site' / 'page.css', docs / 'page.css')
(docs / '.nojekyll').write_text('')
for p in sorted(docs.iterdir()):
    print(f'{p.name:24s} {p.stat().st_size:>10,d} bytes')
