#!/usr/bin/env python3
"""Gaia Archive DataLink probe before Gaia DR4 (run before the archive closes on 2026-11-25).

Asks the live DataLink service for one Gaia DR3 source's epoch photometry under every
DATA_STRUCTURE value (DATAMODEL_STANDARD, INDIVIDUAL, none) and formats, and records what comes
back: HTTP status, size, kind (zip / csv / votable / html / text), zip members and the CSV header.
It also tries Gaia BH3's DR4 epoch astrometry (expected to fail before DR4) and checks whether
gaiadr4 tables answer TAP yet. Standard library only; no analysis.

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
"""
import argparse
import io
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile

BASE = "https://gea.esac.esa.int/data-server/data?"
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
UA = {"User-Agent": "gaiaorb-dr4-probe/0.1 (open research code)"}
DR3_PHOT = 141810650360843008          # G = 10.6, has_epoch_photometry (TAP, 2026-11-22)
BH3_DR3 = 4318465066420528000          # Gaia BH3 (in the DR4 pre-release)


def kind_of(blob):
    if blob[:2] == b"PK":
        return "zip"
    head = blob[:300].lstrip().lower()
    if head.startswith(b"<!doctype html") or head.startswith(b"<html"):
        return "html"
    if head.startswith(b"<?xml") or head.startswith(b"<votable"):
        return "votable"
    first = head.split(b"\n", 1)[0]
    if b"," in first:
        return "csv"
    return "text"


def describe(blob):
    """Summary of one response body (no network)."""
    k = kind_of(blob)
    d = dict(kind=k, bytes=len(blob), head=blob[:300].decode("utf-8", "replace"))
    if k == "zip":
        try:
            with zipfile.ZipFile(io.BytesIO(blob)) as z:
                names = z.namelist()
                d["members"] = names[:20]
                for n in names:
                    if n.lower().endswith(".csv"):
                        d["first_csv_header"] = z.read(n)[:600].decode("utf-8", "replace").split("\n", 1)[0]
                        break
                    if n.lower().endswith((".xml", ".vot", ".ecsv")):
                        d["first_table_head"] = z.read(n)[:600].decode("utf-8", "replace")
                        break
        except zipfile.BadZipFile as e:
            d["error"] = "bad zip: %s" % e
    elif k == "csv":
        d["csv_header"] = blob.split(b"\n", 1)[0][:600].decode("utf-8", "replace")
    return d


def get(url, timeout=180):
    req = urllib.request.Request(url, headers=UA)
    t = time.time()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            blob = r.read()
            out = dict(status=r.status, content_type=r.headers.get("Content-Type"))
    except urllib.error.HTTPError as e:
        blob = e.read() or b""
        out = dict(status=e.code, content_type=e.headers.get("Content-Type") if e.headers else None)
    except Exception as e:
        return dict(status=None, error="%s: %s" % (type(e).__name__, e), seconds=round(time.time() - t, 1))
    out.update(describe(blob))
    out["seconds"] = round(time.time() - t, 1)
    return out


def cases():
    out = []
    for structure in ("DATAMODEL_STANDARD", "INDIVIDUAL", None):
        for fmt in ("CSV", "votable"):
            p = dict(RETRIEVAL_TYPE="EPOCH_PHOTOMETRY", ID="Gaia DR3 %d" % DR3_PHOT, FORMAT=fmt,
                     RELEASE="Gaia DR3")
            if structure:
                p["DATA_STRUCTURE"] = structure
            out.append(("dr3_phot", structure, fmt, BASE + urllib.parse.urlencode(p)))
    for structure in ("DATAMODEL_STANDARD", "INDIVIDUAL"):
        p = dict(RETRIEVAL_TYPE="EPOCH_ASTROMETRY", ID="Gaia DR4 %d" % BH3_DR3, FORMAT="CSV",
                 RELEASE="Gaia DR4", DATA_STRUCTURE=structure)
        out.append(("dr4_bh3_astrometry", structure, "CSV", BASE + urllib.parse.urlencode(p)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/out")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    res = dict(version="0.1.0", started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), datalink=[])
    for what, structure, fmt, url in cases():
        r = get(url)
        r.update(what=what, data_structure=structure, format=fmt, url=url)
        res["datalink"].append(r)
        print(what, structure, fmt, r.get("status"), r.get("kind"), r.get("bytes"), flush=True)
        time.sleep(2)
    for table in ("gaiadr4.gaia_source", "gaiadr3.gaia_source"):
        q = "SELECT TOP 1 source_id FROM %s" % table
        r = get(TAP + "?" + urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q)))
        r.update(what="tap", table=table)
        res.setdefault("tap", []).append(r)
        print("tap", table, r.get("status"), r.get("head", "")[:80].replace("\n", " | "), flush=True)
    with open(os.path.join(args.out_dir, "probe.json"), "w") as f:
        json.dump(res, f, indent=1)
    print("done", flush=True)


if __name__ == "__main__":
    main()
