#!/usr/bin/env python3
"""
fetch_papers.py - weekly literature fetcher for a postdoc search.

Collects the past week's research articles from a fixed set of journals (via PubMed)
plus new bioRxiv preprints, tags each with relevance to the researcher's interests,
and writes JSON that can be pasted into Claude together with weekly_report_prompt.md.

Usage
    python fetch_papers.py                        # last 7 full days, ending yesterday
    python fetch_papers.py --days 14
    python fetch_papers.py --start 2026-09-01 --end 2026-09-07
    python fetch_papers.py --outdir ~/lit-reports

Optional environment variables
    NCBI_API_KEY   free key from your NCBI account; raises the rate limit 3 -> 10 req/s
    NCBI_EMAIL     contact address NCBI asks tools to send (also sent to Crossref)
    CROSSREF_MAILTO  optional separate contact address for Crossref

Requires: Python 3.8+, requests  (pip install requests)
"""

import argparse
import datetime as dt
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET

import requests

# ----------------------------------------------------------------------------- config

# NLM title abbreviation -> display name
JOURNALS = {
    "Nature": "Nature",
    "Nat Genet": "Nature Genetics",
    "Nat Commun": "Nature Communications",
    "Nat Biotechnol": "Nature Biotechnology",
    "Nat Cell Biol": "Nature Cell Biology",
    "Nat Methods": "Nature Methods",
    "Cell": "Cell",
    "Dev Cell": "Developmental Cell",
    "Cell Rep": "Cell Reports",
    "Science": "Science",
    "Sci Adv": "Science Advances",
    "Development": "Development",
}

# PubMed publication types that are not research content
EXCLUDED_PUB_TYPES = {
    "Editorial", "Comment", "News", "Newspaper Article", "Published Erratum",
    "Retraction of Publication", "Interview", "Biography", "Portrait",
    "Introductory Journal Article", "Letter", "Expression of Concern",
}

# bioRxiv categories. Core: kept if they match >= 1 interest theme.
# Extended: kept only if they match >= 2 themes (they are large and mostly off-topic).
BIORXIV_CORE = {"genomics", "developmental biology", "bioinformatics", "systems biology"}
BIORXIV_EXTENDED = {
    "genetics", "cell biology", "molecular biology", "synthetic biology",
    "evolutionary biology", "bioengineering", "neuroscience", "immunology",
    "cancer biology", "plant biology", "zoology",
}

# Researcher interests -> regex patterns (case-insensitive) used for relevance tagging
THEMES = {
    "single_cell_genomics": [
        r"single[- ]cell", r"single[- ]nucle", r"\bsc(?:RNA|ATAC|DNA|Hi-?C)[- ]?seq",
        r"\bsnRNA", r"\bsnATAC", r"multiome", r"spatial(?:ly resolved)? transcriptom",
        r"spatial (?:omics|genomics)", r"cell atlas", r"lineage (?:tracing|recording|barcod)",
        r"perturb-?seq", r"crop-?seq", r"cell[- ]state",
    ],
    "gene_regulatory_networks": [
        r"gene regulatory network", r"\bGRNs?\b", r"transcription factors?",
        r"\benhancers?\b", r"cis-regulatory", r"chromatin accessib", r"\bregulons?\b",
        r"gene regulation", r"transcriptional (?:regulation|control|program)",
        r"\bATAC-seq", r"3D genome", r"chromatin (?:loop|conformation|state)",
        r"pioneer factor",
    ],
    "developmental_biology": [
        r"\bembryo", r"gastrulat", r"organogenesis", r"morphogen", r"\borganoids?\b",
        r"\bgastruloids?", r"\bembryoids?", r"cell fate", r"developmental",
        r"\bsomit", r"neural crest", r"zebrafish", r"drosophila", r"xenopus",
        r"axolotl", r"patterning", r"regenerat", r"lineage specification",
    ],
    "sequence_to_function": [
        r"sequence[- ]to[- ](?:function|expression|activity)", r"deep learning",
        r"neural network", r"language models?", r"foundation models?",
        r"machine learning", r"\bEnformer", r"\bBorzoi", r"AlphaGenome", r"\bBPNet",
        r"ChromBPNet", r"massively parallel reporter", r"\bMPRAs?\b", r"STARR-seq",
        r"variant effect", r"regulatory (?:code|grammar|syntax|logic)",
        r"sequence determinants", r"(?:designed|synthetic) (?:enhancer|promoter)s?",
        r"in silico (?:mutagenesis|design)",
    ],
}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIORXIV_API = "https://api.biorxiv.org/details/biorxiv"
EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works"
USER_AGENT = "weekly-lit-fetch/1.0 (personal literature digest)"

# ----------------------------------------------------------------------------- helpers

_COMPILED = {t: [re.compile(p, re.I) for p in pats] for t, pats in THEMES.items()}
_TAG_RE = re.compile(r"<[^>]+>")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def clean(s):
    if not s:
        return ""
    return " ".join(_TAG_RE.sub(" ", s).split())


def text_of(el):
    return " ".join("".join(el.itertext()).split()) if el is not None else ""


def _ln(tag):
    """Local name of an XML tag (drops any namespace)."""
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else ""


def relevance(text):
    themes, terms = [], set()
    for theme, pats in _COMPILED.items():
        hit = False
        for p in pats:
            m = p.search(text)
            if m:
                hit = True
                terms.add(m.group(0).lower())
        if hit:
            themes.append(theme)
    return {"n_themes": len(themes), "themes": themes, "matched_terms": sorted(terms)}


def http_request(session, method, url, tries=5, **kw):
    kw.setdefault("timeout", 60)
    for i in range(tries):
        try:
            r = session.request(method, url, **kw)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}", response=r)
            r.raise_for_status()
            return r
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as e:
            if i == tries - 1:
                raise
            wait = 2 ** i
            log(f"    retry {i + 1}/{tries - 1} after: {e} (waiting {wait}s)")
            time.sleep(wait)


# ----------------------------------------------------------------------------- PubMed

def _pm_date(article, pubmed_data):
    ad = article.find("ArticleDate")
    try:
        if ad is not None:
            return "%s-%02d-%02d" % (text_of(ad.find("Year")), int(text_of(ad.find("Month"))),
                                     int(text_of(ad.find("Day"))))
        for d in pubmed_data.findall("History/PubMedPubDate"):
            if d.get("PubStatus") == "entrez":
                return "%s-%02d-%02d" % (text_of(d.find("Year")), int(text_of(d.find("Month"))),
                                         int(text_of(d.find("Day"))))
    except (ValueError, TypeError):
        pass
    return ""


def parse_pubmed_xml(xml_bytes):
    root = ET.fromstring(xml_bytes)
    for art in root.findall("PubmedArticle"):
        mc = art.find("MedlineCitation")
        a = mc.find("Article")
        pd = art.find("PubmedData")

        abstract_parts = []
        for at in a.findall("Abstract/AbstractText"):
            label, t = at.get("Label"), text_of(at)
            if t:
                abstract_parts.append(f"{label}: {t}" if label and label.upper() != "UNLABELLED" else t)

        authors = []
        for au in a.findall("AuthorList/Author"):
            coll = au.find("CollectiveName")
            if coll is not None:
                name, individual = text_of(coll), False
            else:
                name = " ".join(x for x in (text_of(au.find("ForeName")), text_of(au.find("LastName"))) if x)
                individual = True
            affs = [text_of(x) for x in au.findall("AffiliationInfo/Affiliation")]
            authors.append({"name": name, "affiliations": affs, "individual": individual})

        individuals = [x for x in authors if x["individual"]] or authors
        last = individuals[-1] if individuals else {"name": "", "affiliations": []}

        doi = ""
        if pd is not None:
            for aid in pd.findall("ArticleIdList/ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = text_of(aid)
        if not doi:
            for e in a.findall("ELocationID"):
                if e.get("EIdType") == "doi":
                    doi = text_of(e)

        pubtypes = [text_of(p) for p in a.findall("PublicationTypeList/PublicationType")]

        yield {
            "pmid": text_of(mc.find("PMID")),
            "journal_abbrev": text_of(mc.find("MedlineJournalInfo/MedlineTA")),
            "title": text_of(a.find("ArticleTitle")),
            "abstract": " ".join(abstract_parts),
            "publication_types": pubtypes,
            "date": _pm_date(a, pd) if pd is not None else "",
            "doi": doi,
            "first_author": individuals[0]["name"] if individuals else "",
            "last_author": {"name": last["name"], "affiliations": last["affiliations"]},
            "n_authors": len(authors),
        }


def fetch_journals(session, start, end, api_key, email, drop_reviews=False):
    base = {"tool": "weekly_lit_fetch"}
    if api_key:
        base["api_key"] = api_key
    if email:
        base["email"] = email
    delay = 0.12 if api_key else 0.4

    jq = " OR ".join(f'"{ta}"[ta]' for ta in JOURNALS)
    term = f'({jq}) AND ("{start:%Y/%m/%d}"[edat] : "{end:%Y/%m/%d}"[edat])'
    log(f"PubMed: searching {len(JOURNALS)} journals, {start} to {end}")
    r = http_request(session, "POST", EUTILS + "esearch.fcgi",
                     data={**base, "db": "pubmed", "term": term, "retmax": "10000", "retmode": "json"})
    ids = r.json()["esearchresult"]["idlist"]
    log(f"  {len(ids)} PubMed records")
    time.sleep(delay)

    stats = {"pubmed_records": len(ids), "excluded_non_research": 0, "excluded_no_abstract": 0,
             "excluded_other_journal": 0, "excluded_reviews": 0, "per_journal": {}}
    records, seen = [], set()
    for i in range(0, len(ids), 200):
        batch = ids[i:i + 200]
        log(f"  fetching records {i + 1}-{i + len(batch)}")
        r = http_request(session, "POST", EUTILS + "efetch.fcgi",
                         data={**base, "db": "pubmed", "id": ",".join(batch), "retmode": "xml"})
        time.sleep(delay)
        for rec in parse_pubmed_xml(r.content):
            if rec["pmid"] in seen:
                continue
            seen.add(rec["pmid"])
            if rec["journal_abbrev"] not in JOURNALS:
                stats["excluded_other_journal"] += 1
                continue
            if EXCLUDED_PUB_TYPES & set(rec["publication_types"]):
                stats["excluded_non_research"] += 1
                continue
            if not rec["abstract"]:
                stats["excluded_no_abstract"] += 1  # mostly News & Views, perspectives, etc.
                continue
            is_review = any("Review" in p for p in rec["publication_types"])
            if is_review and drop_reviews:
                stats["excluded_reviews"] += 1
                continue

            journal = JOURNALS[rec["journal_abbrev"]]
            stats["per_journal"][journal] = stats["per_journal"].get(journal, 0) + 1
            records.append({
                "source": "journal",
                "journal": journal,
                "article_type": "Review" if is_review else "Research article",
                "title": rec["title"],
                "abstract": rec["abstract"],
                "date": rec["date"],
                "doi": rec["doi"],
                "url": f"https://doi.org/{rec['doi']}" if rec["doi"] else f"https://pubmed.ncbi.nlm.nih.gov/{rec['pmid']}/",
                "pmid": rec["pmid"],
                "first_author": rec["first_author"],
                "last_author": rec["last_author"],
                "n_authors": rec["n_authors"],
                "relevance": relevance(rec["title"] + " " + rec["abstract"]),
            })
    log(f"  kept {len(records)} research articles/reviews")
    return records, stats


# ----------------------------------------------------------------------------- bioRxiv

def parse_jats_last_author(xml_bytes):
    """Return {'name', 'affiliations'} for the last author in a bioRxiv JATS XML file."""
    root = ET.fromstring(xml_bytes)
    affs = {}
    for el in root.iter():
        if _ln(el.tag) == "aff":
            label = ""
            for ch in el:
                if _ln(ch.tag) == "label":
                    label = text_of(ch)
            txt = text_of(el)
            if label and txt.startswith(label):
                txt = txt[len(label):].strip(" ,;")
            affs[el.get("id", "")] = txt

    authors = []
    for el in root.iter():
        if _ln(el.tag) == "contrib" and el.get("contrib-type") == "author":
            name, rids = "", []
            for sub in el.iter():
                ln = _ln(sub.tag)
                if ln == "name" and not name:
                    parts = {_ln(x.tag): text_of(x) for x in sub}
                    name = f"{parts.get('given-names', '')} {parts.get('surname', '')}".strip()
                elif ln == "collab" and not name:
                    name = text_of(sub)
                elif ln == "xref" and sub.get("ref-type") == "aff":
                    rids.extend(sub.get("rid", "").split())
            authors.append((name, rids))
    if not authors:
        return None
    name, rids = authors[-1]
    aff_list = [affs[r] for r in rids if r in affs]
    if not aff_list and len(affs) == 1:
        aff_list = list(affs.values())
    return {"name": name, "affiliations": aff_list}


def _surname(name):
    # bioRxiv names look like "Smith, J. A."; Europe PMC names like "Smith JA"
    name = name.strip()
    if "," in name:
        return name.split(",")[0].strip().lower()
    return name.split(" ")[0].strip().lower() if name else ""


def _get_json(session, url, params=None, tries=4):
    """GET and parse JSON; retries when the server returns a non-JSON page (e.g. a bot check)."""
    for i in range(tries):
        r = http_request(session, "GET", url, params=params)
        try:
            return r.json()
        except ValueError:
            snippet = " ".join(r.text[:200].split())
            log(f"    non-JSON response (HTTP {r.status_code}, {r.headers.get('content-type', '?')}): "
                f"{snippet or '<empty body>'}")
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def _biorxiv_candidates(session, start, end):
    raw, cursor, total = [], 0, None
    while True:
        url = f"{BIORXIV_API}/{start:%Y-%m-%d}/{end:%Y-%m-%d}/{cursor}/json"
        data = _get_json(session, url)
        coll = data.get("collection", []) or []
        if total is None:
            try:
                total = int((data.get("messages") or [{}])[0].get("total", 0))
            except (TypeError, ValueError):
                total = 0
            log(f"  {total} bioRxiv records (all versions)")
        raw.extend(coll)
        cursor += len(coll)
        if not coll or cursor >= total:
            break
        time.sleep(0.5)

    out = []
    for p in raw:
        if str(p.get("version", "")) != "1":
            continue
        authors = [a.strip() for a in (p.get("authors") or "").split(";") if a.strip()]
        published = p.get("published")
        out.append({
            "doi": p.get("doi", ""),
            "title": clean(p.get("title")),
            "abstract": clean(p.get("abstract")),
            "date": p.get("date", ""),
            "category": (p.get("category") or "").strip().lower(),
            "authors": authors,
            "last_author_affiliations": [],
            "affiliation_source": None,
            "corresponding_author": p.get("author_corresponding", ""),
            "corresponding_institution": p.get("author_corresponding_institution", ""),
            "published_in": published if published and published != "NA" else None,
            "jats": p.get("jatsxml"),
        })
    return out


def _europepmc_candidates(session, start, end):
    query = (f'SRC:PPR AND PUBLISHER:"bioRxiv" AND '
             f'FIRST_PDATE:[{start:%Y-%m-%d} TO {end:%Y-%m-%d}]')
    params = {"query": query, "resultType": "core", "format": "json",
              "pageSize": "1000", "cursorMark": "*"}
    raw = []
    while True:
        data = _get_json(session, EUROPEPMC, params=params)
        if not raw:
            log(f"  {data.get('hitCount', '?')} bioRxiv preprints in Europe PMC")
        results = (data.get("resultList") or {}).get("result", []) or []
        raw.extend(results)
        nxt = data.get("nextCursorMark")
        if not results or not nxt or nxt == params["cursorMark"]:
            break
        params["cursorMark"] = nxt
        time.sleep(0.3)

    out = []
    for p in raw:
        publisher = ((p.get("bookOrReportDetails") or {}).get("publisher") or "").lower()
        if publisher and "biorxiv" not in publisher:
            continue
        authors, last_affs = [], []
        alist = (p.get("authorList") or {}).get("author", []) or []
        for a in alist:
            authors.append(a.get("fullName") or a.get("collectiveName") or "")
        if alist:
            details = (alist[-1].get("authorAffiliationDetailsList") or {}).get("authorAffiliation", []) or []
            last_affs = [d.get("affiliation", "") for d in details if d.get("affiliation")]
            if not last_affs and alist[-1].get("affiliation"):
                last_affs = [alist[-1]["affiliation"]]
        out.append({
            "doi": p.get("doi", ""),
            "title": clean(p.get("title")),
            "abstract": clean(p.get("abstractText")),
            "date": p.get("firstPublicationDate", ""),
            "category": None,  # Europe PMC doesn't carry bioRxiv subject categories
            "authors": authors,
            "last_author_affiliations": last_affs,
            "affiliation_source": "Europe PMC" if last_affs else None,
            "corresponding_author": "",
            "corresponding_institution": "",
            "published_in": None,
            "jats": None,
        })
    return out


def _jats_abstract(s):
    txt = clean(s)
    return re.sub(r"^(Abstract|Summary)\s*[:.]?\s+", "", txt, flags=re.I)


def _crossref_candidates(session, start, end):
    """bioRxiv preprints via Crossref. bioRxiv registers each preprint's DOI there once, at
    first posting, so filtering on the DOI creation date gives new preprints, not revisions."""
    params = {
        "filter": (f"prefix:10.1101,type:posted-content,"
                   f"from-created-date:{start:%Y-%m-%d},until-created-date:{end:%Y-%m-%d}"),
        "select": "DOI,title,abstract,author,posted,created,group-title,institution,relation",
        "rows": "1000",
        "cursor": "*",
    }
    mailto = os.environ.get("CROSSREF_MAILTO") or os.environ.get("NCBI_EMAIL")
    if mailto:
        params["mailto"] = mailto  # puts requests in Crossref's faster "polite" pool
    items, total = [], None
    while True:
        msg = _get_json(session, CROSSREF, params=params).get("message", {})
        if total is None:
            total = msg.get("total-results", 0)
            log(f"  {total} Crossref records with the bioRxiv/medRxiv DOI prefix")
        batch = msg.get("items", []) or []
        items.extend(batch)
        nxt = msg.get("next-cursor")
        if not batch or not nxt or len(items) >= total:
            break
        params["cursor"] = nxt
        time.sleep(0.3)

    out = []
    for p in items:
        servers = " ".join(i.get("name", "") for i in (p.get("institution") or [])).lower()
        if "medrxiv" in servers:
            continue
        authors, last_affs = [], []
        for a in p.get("author", []) or []:
            name = " ".join(x for x in (a.get("given", ""), a.get("family", "")) if x) or a.get("name", "")
            authors.append(name)
        if p.get("author"):
            last_affs = [x.get("name", "") for x in (p["author"][-1].get("affiliation") or []) if x.get("name")]
        try:
            y, m, d = (p.get("posted") or p.get("created"))["date-parts"][0][:3]
            date = f"{y}-{m:02d}-{d:02d}"
        except (TypeError, KeyError, ValueError, IndexError):
            date = ""
        published = [r.get("id") for r in (p.get("relation") or {}).get("is-preprint-of", []) if r.get("id")]
        group = p.get("group-title")
        out.append({
            "doi": p.get("DOI", ""),
            "title": clean(" ".join(p.get("title") or [])),
            "abstract": _jats_abstract(p.get("abstract", "")),
            "date": date,
            "category": group.strip().lower() if group else None,
            "authors": authors,
            "last_author_affiliations": last_affs,
            "affiliation_source": "Crossref" if last_affs else None,
            "corresponding_author": "",
            "corresponding_institution": "",
            "published_in": published[0] if published else None,
            "jats": None,
        })
    return out


def _europepmc_affiliations(session, dois):
    """Look up last-author affiliations in Europe PMC by DOI, in batches. Returns {doi: [affs]}."""
    found = {}
    for i in range(0, len(dois), 40):
        chunk = dois[i:i + 40]
        query = " OR ".join(f'DOI:"{d}"' for d in chunk)
        try:
            data = _get_json(session, EUROPEPMC, params={"query": query, "resultType": "core",
                                                          "format": "json", "pageSize": "100"}, tries=2)
        except (requests.RequestException, ValueError) as e:
            log(f"    Europe PMC affiliation lookup failed: {e}")
            break
        for p in (data.get("resultList") or {}).get("result", []) or []:
            alist = (p.get("authorList") or {}).get("author", []) or []
            if not alist or not p.get("doi"):
                continue
            det = (alist[-1].get("authorAffiliationDetailsList") or {}).get("authorAffiliation", []) or []
            affs = [d.get("affiliation", "") for d in det if d.get("affiliation")]
            if not affs and alist[-1].get("affiliation"):
                affs = [alist[-1]["affiliation"]]
            if affs:
                found[p["doi"].lower()] = affs
        time.sleep(0.3)
    return found


def fetch_preprints(session, start, end, min_core, min_extended, min_uncategorised, use_jats=True):
    log(f"Preprints: bioRxiv postings {start} to {end}")
    stats = {"source": None, "candidates": 0, "kept": 0, "per_category": {}, "errors": []}
    candidates = None
    for name, fn in (("bioRxiv API", _biorxiv_candidates), ("Crossref", _crossref_candidates),
                     ("Europe PMC", _europepmc_candidates)):
        try:
            log(f"  trying {name}")
            candidates = fn(session, start, end)
            stats["source"] = name
            break
        except (requests.RequestException, ValueError) as e:
            msg = f"{name} failed: {type(e).__name__}: {str(e)[:200]}"
            log("  " + msg)
            stats["errors"].append(msg)
    if candidates is None:
        log("  no preprint source reachable; continuing without preprints")
        return [], stats
    stats["candidates"] = len(candidates)

    kept = {}
    for c in candidates:
        rel = relevance(c["title"] + " " + c["abstract"])
        cat = c["category"]
        if cat is None:
            threshold = min_uncategorised
        elif cat in BIORXIV_CORE:
            threshold = min_core
        elif cat in BIORXIV_EXTENDED:
            threshold = min_extended
        else:
            continue
        if rel["n_themes"] < threshold or not c["doi"]:
            continue
        a = c["authors"]
        kept[c["doi"]] = {
            "source": "preprint",
            "journal": "bioRxiv",
            "category": cat,
            "article_type": "Preprint",
            "title": c["title"],
            "abstract": c["abstract"],
            "date": c["date"],
            "doi": c["doi"],
            "url": f"https://doi.org/{c['doi']}",
            "first_author": a[0] if a else "",
            "last_author": {"name": a[-1] if a else "", "affiliations": c["last_author_affiliations"],
                            "affiliation_source": c["affiliation_source"]},
            "n_authors": len(a),
            "corresponding_author": c["corresponding_author"],
            "corresponding_institution": c["corresponding_institution"],
            "published_in": c["published_in"],
            "relevance": rel,
            "_jats": c["jats"],
        }
        key = cat or "unknown category"
        stats["per_category"][key] = stats["per_category"].get(key, 0) + 1

    records = list(kept.values())
    log(f"  kept {len(records)} relevant preprints (of {len(candidates)} new preprints)")

    # Last-author affiliations, in order of preference for papers still missing one:
    # the source's own data (above), Europe PMC by DOI, the full-text XML on biorxiv.org,
    # and finally the corresponding author's institution when they are the last author.
    missing = [r["doi"] for r in records if not r["last_author"]["affiliations"]]
    if missing:
        log(f"  looking up {len(missing)} last-author affiliations in Europe PMC")
        epmc = _europepmc_affiliations(session, missing)
        for rec in records:
            affs = epmc.get(rec["doi"].lower())
            if affs and not rec["last_author"]["affiliations"]:
                rec["last_author"]["affiliations"] = affs
                rec["last_author"]["affiliation_source"] = "Europe PMC"

    jats_ok, jats_failures = use_jats, 0
    for rec in records:
        jats_url = rec.pop("_jats", None)
        la = rec["last_author"]
        if la["affiliations"]:
            continue
        if jats_ok and jats_url:
            try:
                r = http_request(session, "GET", jats_url, tries=2, timeout=30)
                info = parse_jats_last_author(r.content)
                jats_failures = 0
                if info and info["affiliations"]:
                    la["affiliations"] = info["affiliations"]
                    la["affiliation_source"] = "full-text XML"
            except (requests.RequestException, ET.ParseError):
                jats_failures += 1
                if jats_failures >= 3:
                    log("    full-text XML isn't reachable from here; skipping it")
                    jats_ok = False
            time.sleep(0.3)
        if not la["affiliations"] and rec["corresponding_institution"]:
            sur = _surname(la["name"])
            if sur and sur in rec["corresponding_author"].lower():
                la["affiliations"] = [rec["corresponding_institution"]]
                la["affiliation_source"] = "bioRxiv metadata (last author is corresponding author)"
    for rec in records:
        rec.pop("_jats", None)
    n_aff = sum(1 for r in records if r["last_author"]["affiliations"])
    log(f"  last-author affiliation found for {n_aff}/{len(records)} preprints")
    stats["kept"] = len(records)
    return records, stats


# ----------------------------------------------------------------------------- output

def sort_key(r):
    rel = r["relevance"]
    return (-rel["n_themes"], -len(rel["matched_terms"]), r["journal"], r["title"])


def compact(rec):
    """Smaller record for pasting: drops abstracts of papers that match none of the themes."""
    c = {k: v for k, v in rec.items() if k not in ("pmid", "n_authors")}
    c["relevance"] = {"themes": rec["relevance"]["themes"]}
    c["last_author"] = {"name": rec["last_author"]["name"],
                        "affiliation": (rec["last_author"]["affiliations"] or [""])[0]}
    if rec["relevance"]["n_themes"] == 0:
        c["abstract"] = None
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, default=7, help="window length in days (default 7)")
    ap.add_argument("--start", type=dt.date.fromisoformat, help="start date YYYY-MM-DD")
    ap.add_argument("--end", type=dt.date.fromisoformat, help="end date YYYY-MM-DD (default: yesterday)")
    ap.add_argument("--outdir", default=".", help="where to write the JSON files")
    ap.add_argument("--no-preprints", action="store_true", help="skip bioRxiv")
    ap.add_argument("--no-jats", action="store_true", help="don't fetch full-text XML for preprint affiliations")
    ap.add_argument("--drop-reviews", action="store_true", help="exclude review articles")
    ap.add_argument("--preprint-min-themes-core", type=int, default=1)
    ap.add_argument("--preprint-min-themes-extended", type=int, default=2)
    ap.add_argument("--preprint-min-themes-uncategorised", type=int, default=2,
                    help="threshold when the subject category is unknown (Europe PMC fallback)")
    args = ap.parse_args()

    end = args.end or (dt.date.today() - dt.timedelta(days=1))
    start = args.start or (end - dt.timedelta(days=args.days - 1))
    if start > end:
        ap.error("start date is after end date")

    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT

    articles, jstats = fetch_journals(session, start, end, os.environ.get("NCBI_API_KEY"),
                                      os.environ.get("NCBI_EMAIL"), args.drop_reviews)
    preprints, pstats = [], None
    if not args.no_preprints:
        preprints, pstats = fetch_preprints(session, start, end, args.preprint_min_themes_core,
                                            args.preprint_min_themes_extended,
                                            args.preprint_min_themes_uncategorised, not args.no_jats)
    articles.sort(key=sort_key)
    preprints.sort(key=sort_key)

    out = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "window": {"start": start.isoformat(), "end": end.isoformat()},
        "researcher_interests": list(THEMES),
        "notes": [
            "Journal articles come from PubMed (date = online publication, or date added to PubMed).",
            "Preprints are new bioRxiv postings (from the bioRxiv API, or Crossref / Europe PMC if bioRxiv blocks the request), pre-filtered for relevance to the researcher_interests.",
            "The last author is usually, but not always, the senior/corresponding author.",
            "relevance.themes is a keyword pre-screen, not a judgement; read the abstract.",
        ],
        "stats": {"journals": jstats, "preprints": pstats},
        "journal_articles": articles,
        "preprints": preprints,
    }

    os.makedirs(os.path.expanduser(args.outdir), exist_ok=True)
    stem = os.path.join(os.path.expanduser(args.outdir), f"literature_{start}_{end}")
    full_path, compact_path = stem + ".json", stem + ".compact.json"
    with open(full_path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    small = dict(out, journal_articles=[compact(r) for r in articles],
                 preprints=[compact(r) for r in preprints])
    with open(compact_path, "w", encoding="utf-8") as f:
        json.dump(small, f, ensure_ascii=False, separators=(",", ":"))

    n_rel = sum(1 for r in articles + preprints if r["relevance"]["n_themes"])
    log("")
    log(f"Journal articles: {len(articles)}   Preprints: {len(preprints)}   Matching your themes: {n_rel}")
    for path in (full_path, compact_path):
        size = os.path.getsize(path)
        log(f"  {path}  ({size / 1e6:.2f} MB, roughly {size // 4:,} tokens)")


if __name__ == "__main__":
    main()