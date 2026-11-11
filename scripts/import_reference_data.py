"""Install a local reference-data ZIP after verifying its SHA-256 manifest."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile

ROOT=Path(__file__).resolve().parents[1]

def verified_files(archive):
    with zipfile.ZipFile(archive) as z:
        manifest=json.loads(z.read('manifest.json'))
        if manifest.get('schema')!='nsf-reference-data-v1':raise ValueError('Unknown reference-data schema')
        files=manifest['files']
        if set(z.namelist())!=set(files)|{'manifest.json'}:raise ValueError('Archive entries do not match manifest')
        result={}
        for name,meta in files.items():
            p=PurePosixPath(name)
            if p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name:
                raise ValueError('Unsafe archive path')
            content=z.read(name)
            if len(content)!=meta['bytes'] or hashlib.sha256(content).hexdigest()!=meta['sha256']:
                raise ValueError('Checksum mismatch: '+name)
            result[name]=content
    return result

def install(archive,destination):
    files=verified_files(archive)
    destination=Path(destination).resolve()
    for name,content in files.items():
        p=destination/name
        if p.exists() and (not p.is_file() or p.read_bytes()!=content):
            raise FileExistsError('Different existing file; choose another --dest: '+str(p))
    for name,content in files.items():
        p=destination/name;p.parent.mkdir(parents=True,exist_ok=True)
        if not p.exists():p.write_bytes(content)
    return len(files)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path)
    p.add_argument('--dest',type=Path,default=ROOT/'reference_data')
    a=p.parse_args()
    try:n=install(a.archive,a.dest)
    except (ValueError,KeyError,OSError,zipfile.BadZipFile) as exc:p.error(str(exc))
    print(f'Verified and installed {n} files in {a.dest.resolve()}')

if __name__=='__main__':main()
