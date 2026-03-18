#!/usr/bin/env python3
"""Validate proposed mapping rows against the ES gene mapping schema.

Checks: correct column count, valid gene_id_type, confidence in [0,1],
es_code exists in ontology, no exact duplicates with existing mapping.

Returns structured JSON output for programmatic consumption.

Usage:
    validate_mapping.py \
        --ontology db/ontology/cices_v5.2.tsv \
        --proposed new_rows.tsv \
        [--existing db/mappings/es_gene_mapping.tsv] \
        [--schema db/mappings/schema.json]
"""

import argparse
import csv
import json
import sys

REQUIRED_COLUMNS = ['es_code', 'es_name', 'gene_id', 'gene_id_type',
                    'confidence', 'evidence_source', 'functional_role',
                    'pathway_context', 'notes']

VALID_GENE_ID_TYPES = {'KO', 'COG', 'Pfam', 'EC', 'CAZy', 'custom_hmm'}
VALID_ROLES = {'producer', 'consumer', 'transformer', 'inhibitor'}


def load_ontology_codes(filepath):
    """Load valid ES codes from ontology."""
    codes = set()
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            codes.add(row['es_code'])
    return codes


def load_existing_keys(filepath):
    """Load existing (gene_id, es_code) pairs to check duplicates."""
    keys = set()
    if not filepath:
        return keys
    try:
        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                keys.add((row['gene_id'], row['es_code']))
    except FileNotFoundError:
        pass
    return keys


def validate_row(row, line_num, valid_codes, existing_keys):
    """Validate a single mapping row. Returns list of error strings."""
    errors = []

    # Check required fields are present
    for col in ['es_code', 'gene_id', 'gene_id_type', 'confidence',
                'functional_role', 'evidence_source']:
        if not row.get(col, '').strip():
            errors.append(f'Missing required field: {col}')

    # Validate gene_id_type
    gtype = row.get('gene_id_type', '').strip()
    if gtype and gtype not in VALID_GENE_ID_TYPES:
        errors.append(f'Invalid gene_id_type: {gtype} (allowed: {", ".join(sorted(VALID_GENE_ID_TYPES))})')

    # Validate confidence
    try:
        conf = float(row.get('confidence', ''))
        if conf < 0 or conf > 1:
            errors.append(f'Confidence {conf} out of range [0, 1]')
    except (ValueError, TypeError):
        errors.append(f'Invalid confidence value: {row.get("confidence", "")}')

    # Validate functional_role
    role = row.get('functional_role', '').strip()
    if role and role not in VALID_ROLES:
        errors.append(f'Invalid functional_role: {role} (allowed: {", ".join(sorted(VALID_ROLES))})')

    # Validate es_code against ontology
    es_code = row.get('es_code', '').strip()
    if es_code and es_code not in valid_codes:
        errors.append(f'ES code {es_code} not found in ontology')

    # Check for exact duplicate
    key = (row.get('gene_id', '').strip(), es_code)
    if key in existing_keys:
        errors.append(f'Duplicate: gene_id={key[0]} + es_code={key[1]} already exists')

    return errors


def main():
    parser = argparse.ArgumentParser(description='Validate proposed ES mappings')
    parser.add_argument('--ontology', required=True, help='ES ontology TSV')
    parser.add_argument('--proposed', required=True, help='Proposed mappings TSV')
    parser.add_argument('--existing', help='Existing mapping table (for dedup check)')
    parser.add_argument('--schema', help='Schema JSON (for reference)')
    parser.add_argument('--json-output', help='Write structured JSON results to file')
    args = parser.parse_args()

    valid_codes = load_ontology_codes(args.ontology)
    existing_keys = load_existing_keys(args.existing)

    results = []
    total = 0
    passed = 0

    with open(args.proposed) as f:
        reader = csv.DictReader(f, delimiter='\t')

        # Check columns
        if reader.fieldnames:
            missing = set(REQUIRED_COLUMNS) - set(reader.fieldnames)
            if missing:
                print(f'WARNING: Missing columns in proposed file: {missing}',
                      file=sys.stderr)

        for i, row in enumerate(reader, start=2):  # line 1 is header
            total += 1
            errors = validate_row(row, i, valid_codes, existing_keys)
            status = 'pass' if not errors else 'fail'
            if not errors:
                passed += 1

            results.append({
                'line': i,
                'gene_id': row.get('gene_id', ''),
                'es_code': row.get('es_code', ''),
                'status': status,
                'errors': errors,
            })

    # Output
    summary = {
        'total': total,
        'passed': passed,
        'failed': total - passed,
        'results': results,
    }

    if args.json_output:
        with open(args.json_output, 'w') as f:
            json.dump(summary, f, indent=2)

    # Also print human-readable summary
    print(f'Validation: {passed}/{total} rows passed', file=sys.stderr)
    for r in results:
        if r['status'] == 'fail':
            print(f"  Line {r['line']}: {r['gene_id']} → {r['es_code']}: "
                  f"{'; '.join(r['errors'])}", file=sys.stderr)

    # Exit with error code if any failures
    sys.exit(0 if passed == total else 1)


if __name__ == '__main__':
    main()
