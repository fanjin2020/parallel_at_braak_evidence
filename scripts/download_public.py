"""Fetch public source data explicitly; the analysis launcher never calls this module."""
import argparse
from datetime import datetime, timezone
import gzip
import http.client
import json
import os
from pathlib import Path
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import zlib
from data_sources import SOURCES, COUNTS, matrix_sample_ids, cel_selection, gaf_release, load_cached_public
from source_io import sha256

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_HOSTS = frozenset({'ftp.ncbi.nlm.nih.gov', 'ftp.ebi.ac.uk', 'media.springernature.com',
    'static-content.springer.com', 'release.geneontology.org', 'www.ebi.ac.uk'})


def public_url(url):
    p = urllib.parse.urlsplit(url)
    if p.scheme != 'https' or p.hostname not in PUBLIC_HOSTS or p.username or p.password or p.port not in (None, 443):
        raise ValueError('Not an allowed public-data HTTPS source: '+url)
    return url


class PublicRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, newurl):
        public_url(newurl)
        return super().redirect_request(request, response, code, message, headers, newurl)


def open_public(url, headers=None):
    if os.environ.get('JAD_OFFLINE') == '1':
        raise RuntimeError('Offline analysis cannot download. Run download.sh separately.')
    request = urllib.request.Request(public_url(url), headers={
        'User-Agent': 'parallel-at-braak-evidence/5.0', **(headers or {})})
    return urllib.request.build_opener(PublicRedirect()).open(request, timeout=120)


def atomic_json(path, value):
    temp = path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temp.replace(path)


def file_format(path, name=None):
    name = name or path.name
    if not path.stat().st_size:
        raise ValueError('Empty source file: '+name)
    if name.endswith('.gz'):
        with gzip.open(path, 'rb') as stream:
            while stream.read(1024*1024):
                pass
    elif name.endswith('.xls'):
        with path.open('rb') as stream:
            if stream.read(8) != bytes.fromhex('d0cf11e0a1b11ae1'):
                raise ValueError('Expected XLS, not an HTML error response: '+name)
    elif name.endswith('.xlsx'):
        with zipfile.ZipFile(path) as archive:
            if archive.testzip() is not None or 'xl/workbook.xml' not in archive.namelist():
                raise ValueError('Invalid XLSX: '+name)
    elif name.endswith('.h5'):
        with path.open('rb') as stream:
            if stream.read(8) != b'\x89HDF\r\n\x1a\n':
                raise ValueError('Invalid HDF5 source: '+name)


def local_identity(path, expected=None):
    expected = expected or {}
    digest = sha256(path)
    if expected.get('sha256') and digest != expected['sha256']:
        raise ValueError('Existing snapshot differs: '+str(path)+'. It will not be overwritten.')
    if expected.get('bytes') is not None and path.stat().st_size != int(expected['bytes']):
        raise ValueError('Source size differs: '+str(path))
    return dict(bytes=path.stat().st_size, sha256=digest)


def fetch(url, path, expected=None, platform_table=False):
    """Reuse complete files; resume partial downloads; never replace an existing snapshot."""
    public_url(url)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path = path.with_name(path.name+'.receipt.json')
    if expected is None and receipt_path.is_file():
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        if receipt.get('url') != url:
            raise ValueError('Cached source URL changed: '+str(path))
        expected = receipt
    if path.exists():
        if not path.is_file():
            raise ValueError('Expected a source file: '+str(path))
        identity = local_identity(path, expected)
        file_format(path)
        print('Using existing public file: '+str(path), flush=True)
        return path, dict(**identity, origin='existing_local_snapshot')
    part = path.with_name(path.name+'.part')
    for attempt in range(3):
        try:
            if platform_table:
                with open_public(url) as response, gzip.GzipFile(fileobj=response) as reader, part.open('wb') as raw:
                    with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as writer:
                        for line in reader:
                            writer.write(line)
                            if line.startswith(b'!platform_table_end'):
                                break
                        else:
                            raise ValueError('SOFT platform table is incomplete')
            else:
                offset = part.stat().st_size if part.is_file() else 0
                with open_public(url, {'Range': f'bytes={offset}-'} if offset else {}) as response:
                    append = bool(offset and response.status == 206)
                    total = response.headers.get('Content-Length')
                    total = int(total)+(offset if append else 0) if total else None
                    if response.status == 206:
                        match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                        if not match or int(match[1]) != (offset if append else 0):
                            raise ValueError('Unexpected HTTP partial-response range')
                        total = int(match[3])
                    last = time.monotonic()
                    with part.open('ab' if append else 'wb') as stream:
                        while block := response.read(1024*1024):
                            stream.write(block)
                            if time.monotonic()-last >= 15:
                                print(f'{path.name}: {stream.tell():,}/{total or "?"} bytes', flush=True)
                                last = time.monotonic()
                if total is not None and part.stat().st_size != total:
                    raise IOError('Incomplete download: '+path.name)
            file_format(part, path.name)
            identity = local_identity(part, expected)
            part.replace(path)
            receipt = dict(url=url, **identity, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                           kind='upstream_platform_table_projection' if platform_table else 'complete_upstream_file')
            atomic_json(receipt_path, receipt)
            return path, dict(**identity, origin='downloaded_public_source')
        except (OSError, ValueError, EOFError, http.client.IncompleteRead, zlib.error) as error:
            if platform_table or isinstance(error, (ValueError, EOFError, zlib.error)) or (
                    isinstance(error, urllib.error.HTTPError) and error.code == 416):
                part.unlink(missing_ok=True)
            if attempt == 2:
                raise
            print(f'Retrying {path.name}: {error}', flush=True)
            time.sleep(2**attempt)


def existing_path(cache, relative, saved=None):
    target = cache/relative
    if target.is_file():
        return target
    if saved and Path(saved).is_file():
        return Path(saved)
    return target


def fetch_main(cache, overrides):
    cache = cache.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    manifest = cache/'public_sources_manifest.json'
    prior = json.loads(manifest.read_text(encoding='utf-8')) if manifest.is_file() else {}
    if prior and 'paths' not in prior:
        raise ValueError('This directory already contains a different public-data manifest; choose a separate main cache.')
    records = {r['key']: r for r in prior.get('sources', [])}
    saved_paths = prior.get('paths', {})
    paths, result = {}, []
    for key, relative, url in SOURCES:
        old = records.get(key)
        if overrides.get(key):
            path = overrides[key].resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            identity = local_identity(path, old)
            file_format(path)
            identity['origin'] = 'user_supplied_public_snapshot'
        else:
            path = existing_path(cache, relative, saved_paths.get(key))
            path, identity = fetch((old or {}).get('source_url') or url, path, old,
                                   platform_table=key == 'external_AD_platform')
        paths[key] = str(path.resolve())
        result.append(dict(key=key, path=paths[key], source_url=(old or {}).get('source_url') or url,
                           kind='upstream_platform_table_projection' if key == 'external_AD_platform' else 'complete_upstream_file',
                           **identity))
    release = gaf_release(paths['go_gaf'])
    ontology_url = f'https://release.geneontology.org/{release}/ontology/go-basic.obo'
    old = records.get('go_ontology')
    ontology = overrides.get('go_ontology')
    if ontology is not None:
        ontology = ontology.resolve()
        identity = dict(**local_identity(ontology, old), origin='user_supplied_public_snapshot')
    else:
        ontology = existing_path(cache, 'annotations/go-basic.obo', saved_paths.get('go_ontology'))
        if old and not ontology.is_file() and old.get('kind') == 'version_checked_public_GO_graph_projection':
            raise FileNotFoundError('Restore the recorded GO projection or use a new cache with matching --go-gaf and --go-ontology.')
        ontology, identity = fetch(ontology_url, ontology, old)
    versions = [line.strip() for line in ontology.read_text(encoding='utf-8').splitlines()
                if line.startswith('data-version:')]
    if len(versions) != 1 or release not in versions[0]:
        raise ValueError('The GO ontology must match the release declared by the GAF.')
    paths['go_ontology'] = str(ontology.resolve())
    result.append(dict(key='go_ontology', path=paths['go_ontology'],
                       source_url=(old or {}).get('source_url') or ontology_url,
                       expected_matching_release_url=ontology_url,
                       kind=(old or {}).get('kind', 'complete_upstream_file'), **identity))
    ids = matrix_sample_ids(paths['discovery_matrix'])
    selected = cel_selection(ids, paths['cel_file_list'])
    cel_dir = cache/'GSE131617/CEL'
    if not cel_dir.is_dir() and Path(saved_paths.get('cel_dir', cel_dir)).is_dir():
        cel_dir = Path(saved_paths['cel_dir'])
    old_cels = {r['sample_id']: r for r in prior.get('cel_files', [])}
    cels = []
    for index, row in enumerate(selected, 1):
        name, sample = row['Name'], row['Name'].split('_')[0]
        url = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{sample[:-3]}nnn/{sample}/suppl/{name}'
        expected = dict(old_cels.get(sample, {}), bytes=int(row['Size']))
        path, identity = fetch(url, cel_dir/name, expected)
        cels.append(dict(sample_id=sample, filename=str(path.resolve()), source_url=url, **identity))
        print(f'CEL inputs: {index}/{len(selected)}', flush=True)
    paths['cel_dir'] = str(cel_dir.resolve())
    record = dict(schema_version=1, restricted_data_downloaded=False, paths=paths, sources=result,
                  discovery_samples=len(ids), cel_download_requested=True, cel_files=cels)
    atomic_json(manifest, record)
    load_cached_public(cache, progress=False)
    print('Public input manifest: '+str(manifest), flush=True)


def fetch_single_nucleus(cache):
    cache = cache.resolve()
    manifest=cache/'public_sources_manifest.json'
    if manifest.is_file():
        prior=json.loads(manifest.read_text(encoding='utf-8'))
        if 'paths' in prior or any('key' in r for r in prior.get('sources',[])):
            raise ValueError('This directory contains the main public-data manifest; choose a separate single-nucleus cache.')
    records = []
    for name, digest in COUNTS.items():
        url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/'+name
        path, identity = fetch(url, cache/name, dict(sha256=digest))
        records.append(dict(filename=name, path=str(path.resolve()), url=url, **identity))
    atomic_json(cache/'public_sources_manifest.json', dict(sources=records, restricted_data_downloaded=False))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache-dir', type=Path, default=ROOT/'cache/public')
    p.add_argument('--public-cache', type=Path, help='GSE174367 cache; defaults to <cache-dir>/GSE174367')
    p.add_argument('--only', choices=['all', 'main', 'single-nucleus'], default='all')
    p.add_argument('--go-gaf', type=Path)
    p.add_argument('--go-ontology', type=Path)
    p.add_argument('--gene-info', type=Path)
    a = p.parse_args()
    main_cache=a.cache_dir.resolve()
    single_cache=(a.public_cache or a.cache_dir/'GSE174367').resolve()
    if main_cache==single_cache or (main_cache.exists() and single_cache.exists() and os.path.samefile(main_cache,single_cache)):
        raise ValueError('--cache-dir and --public-cache must be different directories. No download was started.')
    if a.only in ('all', 'main'):
        fetch_main(a.cache_dir, dict(go_gaf=a.go_gaf, go_ontology=a.go_ontology, gene_info=a.gene_info))
    if a.only in ('all', 'single-nucleus'):
        fetch_single_nucleus(single_cache)


if __name__ == '__main__':
    main()

