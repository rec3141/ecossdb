#!/usr/bin/env python3
"""Normalize annotation inputs into a unified format.

Accepts multiple input formats and produces a standardized TSV with columns:
protein_id, contig_id, KO, COG, EC, Pfam, CAZy, description

Supported formats:
  - danaseq merged_annotations.tsv (auto-detected)
  - kofamscan detail output
  - eggNOG-mapper annotations
  - GFF3 with functional annotations
  - Simple KO list (one KO per line, or protein_id<TAB>KO)
  - dbCAN overview output

Usage:
    normalize_annotations.py --input annotations.tsv --format auto --output normalized.tsv
"""

import argparse
import csv
import os
import re
import sys


def detect_format(filepath):
    """Auto-detect annotation file format from header/content."""
    with open(filepath) as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip() if first_line else ''

    header_lower = first_line.lower()

    # danaseq merged format
    if 'protein_id' in header_lower and 'ko' in header_lower and 'contig_id' in header_lower:
        return 'danaseq'

    # eggNOG-mapper
    if first_line.startswith('#') and 'eggNOG' in first_line:
        return 'eggnog'
    if 'seed_ortholog' in header_lower or 'eggnog' in header_lower.replace('-', ''):
        return 'eggnog'

    # kofamscan
    if header_lower.startswith('#') or (second_line and re.match(r'^\*?\s+K\d{5}', second_line)):
        return 'kofamscan'
    if re.match(r'^\*?\s+K\d{5}', first_line):
        return 'kofamscan'

    # GFF
    if '##gff-version' in first_line or first_line.count('\t') >= 7:
        # Check if it looks like GFF (9 columns, col 3 is feature type)
        parts = first_line.split('\t')
        if len(parts) >= 9 and parts[2] in ('CDS', 'gene', 'mRNA'):
            return 'gff'

    # dbCAN
    if 'hmm' in header_lower and ('cazy' in header_lower or 'subfamily' in header_lower):
        return 'dbcan'

    # Simple KO list
    if re.match(r'^K\d{5}', first_line):
        return 'ko_list'
    parts = first_line.split('\t')
    if len(parts) == 2 and re.match(r'^K\d{5}', parts[1]):
        return 'ko_list'

    return 'danaseq'  # fallback


def clean_ko(value):
    """Strip ko: prefix and clean KO identifiers."""
    if not value:
        return ''
    # Handle comma-separated values with ko: prefix
    kos = []
    for v in value.split(','):
        v = v.strip().replace('ko:', '')
        if v:
            kos.append(v)
    return ','.join(kos)


def parse_danaseq(filepath):
    """Parse danaseq merged_annotations.tsv format."""
    records = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        cols = set(reader.fieldnames or [])
        for row in reader:
            records.append({
                'protein_id': row.get('protein_id', ''),
                'contig_id': row.get('contig_id', ''),
                'KO': clean_ko(row.get('KO', row.get('ko', ''))),
                'COG': row.get('COG', row.get('COG_category', row.get('cog', ''))),
                'EC': row.get('EC', row.get('ec', '')),
                'Pfam': row.get('Pfam', row.get('PFAMs', row.get('pfam', ''))),
                'CAZy': row.get('CAZy', row.get('cazy', row.get('CAZy_family', ''))),
                'description': row.get('description', row.get('product', '')),
            })
    return records


def parse_kofamscan(filepath):
    """Parse kofamscan detail/tabular output."""
    records = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # kofamscan format: significance gene_name KO threshold score E-value KO_definition
            parts = line.split()
            if len(parts) < 5:
                continue
            sig = parts[0]
            protein_id = parts[1] if sig in ('*', '') else parts[0]
            ko = parts[2] if sig in ('*', '') else parts[1]
            if not ko.startswith('K'):
                continue
            # Only keep significant hits
            if sig != '*' and len(parts) > 2:
                # Try tab-separated format
                tparts = line.split('\t')
                if len(tparts) >= 3:
                    protein_id = tparts[0].strip().lstrip('* ')
                    ko = tparts[1].strip()
                    if not ko.startswith('K'):
                        continue
            contig_id = '_'.join(protein_id.rsplit('_', 1)[:-1]) if '_' in protein_id else protein_id
            records.append({
                'protein_id': protein_id,
                'contig_id': contig_id,
                'KO': ko,
                'COG': '', 'EC': '', 'Pfam': '', 'CAZy': '',
                'description': ' '.join(parts[6:]) if len(parts) > 6 else '',
            })
    return records


def parse_eggnog(filepath):
    """Parse eggNOG-mapper annotations."""
    records = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                # Check for header
                if line.startswith('#query'):
                    header = line.lstrip('#').strip().split('\t')
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            protein_id = parts[0]
            contig_id = '_'.join(protein_id.rsplit('_', 1)[:-1]) if '_' in protein_id else protein_id

            # eggNOG v2 format columns vary, find by position or name
            ko = ''
            cog = ''
            ec = ''
            pfam = ''
            desc = ''

            # Standard eggNOG-mapper v2 output positions
            if len(parts) > 11:
                ko = parts[11] if len(parts) > 11 else ''  # KEGG_ko
                ec = parts[10] if len(parts) > 10 else ''   # EC
                cog = parts[6] if len(parts) > 6 else ''    # COG_category
                pfam = parts[20] if len(parts) > 20 else '' # PFAMs
                desc = parts[7] if len(parts) > 7 else ''   # Description

            # Clean KO field (may have 'ko:' prefix)
            ko = ','.join(k.replace('ko:', '') for k in ko.split(',') if k.startswith(('ko:', 'K')))

            records.append({
                'protein_id': protein_id,
                'contig_id': contig_id,
                'KO': ko, 'COG': cog, 'EC': ec, 'Pfam': pfam, 'CAZy': '',
                'description': desc,
            })
    return records


def parse_gff(filepath):
    """Parse GFF3 with functional annotations in attributes."""
    records = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CDS':
                continue

            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attrs[k.strip()] = v.strip()

            protein_id = attrs.get('ID', attrs.get('locus_tag', ''))
            contig_id = parts[0]

            records.append({
                'protein_id': protein_id,
                'contig_id': contig_id,
                'KO': attrs.get('KO', attrs.get('ko', '')),
                'COG': attrs.get('COG', attrs.get('cog', '')),
                'EC': attrs.get('eC_number', attrs.get('EC', attrs.get('ec', ''))),
                'Pfam': attrs.get('Pfam', attrs.get('pfam', '')),
                'CAZy': attrs.get('CAZy', attrs.get('cazy', '')),
                'description': attrs.get('product', attrs.get('Name', '')),
            })
    return records


def parse_ko_list(filepath):
    """Parse simple KO list (one per line or protein_id<TAB>KO)."""
    records = []
    counter = 0
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                protein_id = parts[0]
                ko = parts[1]
            else:
                counter += 1
                protein_id = f'protein_{counter:06d}'
                ko = parts[0]

            if not re.match(r'^K\d{5}', ko):
                continue

            records.append({
                'protein_id': protein_id,
                'contig_id': '',
                'KO': ko,
                'COG': '', 'EC': '', 'Pfam': '', 'CAZy': '',
                'description': '',
            })
    return records


def parse_dbcan(filepath):
    """Parse dbCAN overview output."""
    records = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            protein_id = row.get('Gene ID', row.get('#ofTools', ''))
            # dbCAN overview has HMMER, DIAMOND, and Hotpep columns
            cazy_hits = []
            for tool in ['HMMER', 'DIAMOND', 'dbCAN_sub']:
                val = row.get(tool, '').strip()
                if val and val != '-':
                    cazy_hits.append(val.split('(')[0])

            if not cazy_hits:
                continue

            contig_id = '_'.join(protein_id.rsplit('_', 1)[:-1]) if '_' in protein_id else protein_id
            records.append({
                'protein_id': protein_id,
                'contig_id': contig_id,
                'KO': '', 'COG': '', 'EC': '', 'Pfam': '',
                'CAZy': ','.join(sorted(set(cazy_hits))),
                'description': '',
            })
    return records


PARSERS = {
    'danaseq': parse_danaseq,
    'kofamscan': parse_kofamscan,
    'eggnog': parse_eggnog,
    'gff': parse_gff,
    'ko_list': parse_ko_list,
    'dbcan': parse_dbcan,
}


def main():
    parser = argparse.ArgumentParser(description='Normalize annotation inputs')
    parser.add_argument('--input', required=True, nargs='+', help='Input file(s)')
    parser.add_argument('--format', default='auto',
                        choices=['auto'] + list(PARSERS.keys()),
                        help='Input format (default: auto-detect)')
    parser.add_argument('--output', required=True, help='Output normalized TSV')
    args = parser.parse_args()

    all_records = []
    for filepath in args.input:
        fmt = args.format
        if fmt == 'auto':
            fmt = detect_format(filepath)
            print(f'Auto-detected format: {fmt} for {os.path.basename(filepath)}',
                  file=sys.stderr)

        parse_fn = PARSERS[fmt]
        records = parse_fn(filepath)
        all_records.extend(records)
        print(f'Parsed {len(records)} records from {os.path.basename(filepath)}',
              file=sys.stderr)

    # Deduplicate by protein_id (keep first occurrence, merge annotations)
    merged = {}
    for rec in all_records:
        pid = rec['protein_id']
        if pid in merged:
            existing = merged[pid]
            for col in ('KO', 'COG', 'EC', 'Pfam', 'CAZy'):
                if rec[col] and not existing[col]:
                    existing[col] = rec[col]
            if rec['description'] and not existing['description']:
                existing['description'] = rec['description']
        else:
            merged[pid] = rec

    cols = ['protein_id', 'contig_id', 'KO', 'COG', 'EC', 'Pfam', 'CAZy', 'description']
    with open(args.output, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for rec in merged.values():
            f.write('\t'.join(str(rec.get(c, '')) for c in cols) + '\n')

    print(f'Total: {len(merged)} normalized annotations', file=sys.stderr)


if __name__ == '__main__':
    main()
