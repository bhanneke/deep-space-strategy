"""Prose gate for the site pages: the FIRM deck's slopscan (Security_Contagion_FIRM/Presentations/slopscan.py)
applied to HTML pages, plus a passive-voice check. Usage: python3 site/slopscan.py site/index.html [...]
Exit 1 on any blocking hit. The regexes are a floor; read every sentence as well."""
import re, pathlib, html, sys

PATTERNS = [
 ("not X, but Y",      r'\bnot\s+[^.,;]{1,45},\s*but\b'),
 ("bullet opens on a demonstrative", r'<li[^>]*>\s*(?:<b>|<strong>)?\s*(?:That|This|These|Those|It|They)\b'),
 ("neither X nor Y",    r'\bneither\b[^.<;]{2,50}\bnor\b'),
 ("It is not X. It is Y", r'\b(is|are|was|were)\s+not\s+[^.;<]{2,45}\.\s*(It|They)\s+(is|are|was|were)\b'),
 ("appositive X, not Y", r'[\w\)\}"][,;]\s+not\s+[^.<;]{2,60}'),
 ("dash - not Y",       r'—\s*not\s+[^.<;]{2,60}'),
 ("X and not Y",        r'\band not\b[^.<;]{2,60}'),
 ("rather than",        r'\brather than\b'),
 ("instead of",         r'\binstead of\b'),
 ("X is not Y; it is Z", r'\b(is|are)\s+not\s+[^.;]{1,45};\s*(it|they)\s+(is|are)\b'),
 ("doesn't X, it Ys",  r"\b(does|do|did|is|are)n[’']t\s+[^.,;]{1,45},\s*(it|they)\b"),
 ("Against this background", r'\b[Aa]gainst this back(ground|drop)\b'),
 ("rule-of-three A, B, and C", r',\s[^,.;:<]{3,40},\s+and\s+[^,.;:<]{3,40}'),
 ("colon + 3-item list", r':\s+[a-z][^.;<]{0,60},[^.;<]{0,60},[^.;<]{0,60}\.'),
 ("mirrored can X; cannot Y", r'\bcan\s+[^.;]{1,40};\s*(it\s+)?cannot\b'),
 ("slogan pair No X. No Y.", r'\bNo\s+[A-Za-z][^.]{0,30}\.\s+No\s+[A-Za-z]'),
 ("banned vocab", r'\b(battle-tested|first principles|right-sized|data-driven|research-backed|simulation-first|bespoke|actually|landscape|game-changer|cutting-edge|groundbreaking|delve)\b'),
 ("rhetorical question", r'<h[12][^>]*>[^<]*\?'),
 ("vague referent (all/both of them)", r'\b(all|both|each|none) of them\b'),
 ("bare demonstrative object", r'\b(hides|shows|explains|means|does|reads|call|calls) this\b'),
 ("and no X / and not Y", r'\b(is|was|are|were)\s+[^.;<]{2,50}\band (no|not)\b'),
 ("produces/changes nothing", r'\b(produces|changes|moves|conveys) nothing\b'),
 ("swipe verbs", r'\b(beats?|trumps?|outperforms?)\b'),
 ("most/overwhelmingly without a count", r'\b(overwhelmingly|the vast majority)\b'),
 ("sentence opens on a demonstrative", r'(?:^|[.>]\s+)(?:That|This|These|Those)\s+(?:makes|means|shows|is|are)\b'),
 ("metaphor: abstract noun doing a physical verb", r'\b(redraws?|reorders?|carries|binds|bends|sets the terms|holds the|strongholds?)\b'),
 ("passive voice", r'\b(?:is|are|was|were|be|been|being)\s+(?:\w+ly\s+)?(?:\w+ed|built|made|done|run|set|given|taken|shown|seen|known|found|held|led|drawn|chosen)\b'),
]
HARD = {name for name, _ in PATTERNS} - {"rule-of-three A, B, and C", "colon + 3-item list", "metaphor: abstract noun doing a physical verb"}

total_hard = 0
for f in sys.argv[1:]:
    s = pathlib.Path(f).read_text()
    body = re.sub(r'<script\b.*?</script>|<style\b.*?</style>|<title>.*?</title>|<meta[^>]*>|\$\$.*?\$\$', '', s, flags=re.S)
    body = html.unescape(body)
    print(f"=== {f}")
    for t in re.findall(r'<h[12][^>]*>(.*?)</h[12]>', body, re.S):
        t = re.sub(r'<[^>]+>', '', t).strip()
        flag = (" <-- TWO SENTENCES" if re.search(r'\.\s+\S', t) else "") + (" <-- COLON" if ':' in t else "")
        if flag:
            print("  heading:", t, flag)
    for name, pat in PATTERNS:
        hits = [m.group(0).replace("\n", " ")[:110] for m in re.finditer(pat, body)]
        if hits:
            hard = name in HARD
            total_hard += len(hits) if hard else 0
            print(f"  [{name}]{'  <<< BLOCKS' if hard else '  (review)'} {len(hits)}")
            for h in hits:
                print("      ", h)
print(f"\nblocking hits: {total_hard}")
sys.exit(1 if total_hard else 0)
