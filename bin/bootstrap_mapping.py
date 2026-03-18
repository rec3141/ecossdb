#!/usr/bin/env python3
"""Bootstrap the ES gene mapping table from KEGG-Decoder modules and FOAM.

Generates the initial es_gene_mapping.tsv by:
1. Mapping KEGG module definitions to CICES Regulation & Maintenance codes
2. Mapping FOAM L1 functional categories to CICES codes
3. Cross-referencing to boost confidence where both sources agree

Usage:
    bootstrap_mapping.py \
        --foam /data/scratch/refdbs/FOAM/FOAM-onto_rel1.tsv \
        --ontology db/ontology/cices_v5.2.tsv \
        --output db/mappings/es_gene_mapping.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# KEGG module → CICES mapping
# Maps biogeochemical KEGG modules to ecosystem service codes and roles.
# Based on KEGG-Decoder categories from danaseq kegg_module_completeness.py.
# ---------------------------------------------------------------------------

KEGG_MODULE_ES = {
    # Nitrogen metabolism → Water quality regulation
    'M00175': {  # Nitrogen fixation
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'producer',
        'confidence': 0.85,
        'pathway_context': 'Nitrogen fixation (N2→NH4)',
    },
    'M00528': {  # Nitrification
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'transformer',
        'confidence': 0.80,
        'pathway_context': 'Nitrification (NH3→NO3)',
    },
    'M00529': {  # Denitrification
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'producer',
        'confidence': 0.85,
        'pathway_context': 'Denitrification (NO3→N2)',
    },
    'M00530': {  # DNRA
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'transformer',
        'confidence': 0.75,
        'pathway_context': 'DNRA (NO3→NH4)',
    },
    'M00804': {  # Anammox
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'producer',
        'confidence': 0.85,
        'pathway_context': 'Anammox (NH4+NO2→N2)',
    },
    'M00531': {  # Assimilatory nitrate reduction
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'consumer',
        'confidence': 0.65,
        'pathway_context': 'Assimilatory nitrate reduction (NO3→NH4)',
    },

    # Sulfur metabolism → Water quality / Waste mediation
    'M00595': {  # Sox thiosulfate oxidation
        'es_code': '2.1.1.1',
        'es_name': 'Decomposing wastes or polluting substances',
        'functional_role': 'transformer',
        'confidence': 0.75,
        'pathway_context': 'Thiosulfate oxidation (Sox pathway)',
    },
    'M00596': {  # Dissimilatory sulfate reduction
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'transformer',
        'confidence': 0.75,
        'pathway_context': 'Dissimilatory sulfate reduction (Dsr)',
    },
    'M00176': {  # Assimilatory sulfate reduction
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'consumer',
        'confidence': 0.60,
        'pathway_context': 'Assimilatory sulfate reduction',
    },

    # Methane metabolism → Climate regulation
    'M00567': {  # Methanogenesis
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'inhibitor',
        'confidence': 0.90,
        'pathway_context': 'Methanogenesis (CO2→CH4)',
    },
    'M00174': {  # Methane oxidation
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.90,
        'pathway_context': 'Methane oxidation (CH4→CO2)',
    },
    'M00346': {  # Formaldehyde assimilation
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.70,
        'pathway_context': 'Formaldehyde assimilation (methylotrophy)',
    },

    # Carbon fixation → Climate regulation
    'M00165': {  # Calvin cycle
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.85,
        'pathway_context': 'Calvin cycle (CO2 fixation)',
    },
    'M00173': {  # rTCA
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.80,
        'pathway_context': 'Reductive TCA cycle (CO2 fixation)',
    },
    'M00377': {  # Wood-Ljungdahl
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.80,
        'pathway_context': 'Wood-Ljungdahl pathway (CO2 fixation)',
    },
    'M00376': {  # 3-HP
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.75,
        'pathway_context': '3-Hydroxypropionate bi-cycle (CO2 fixation)',
    },

    # Xenobiotics degradation → Waste decomposition
    'M00548': {  # Benzene degradation
        'es_code': '2.1.1.1',
        'es_name': 'Decomposing wastes or polluting substances',
        'functional_role': 'producer',
        'confidence': 0.85,
        'pathway_context': 'Benzene degradation',
    },
    'M00568': {  # Catechol ortho-cleavage
        'es_code': '2.1.1.1',
        'es_name': 'Decomposing wastes or polluting substances',
        'functional_role': 'producer',
        'confidence': 0.80,
        'pathway_context': 'Catechol ortho-cleavage (xenobiotic degradation)',
    },
    'M00569': {  # Catechol meta-cleavage
        'es_code': '2.1.1.1',
        'es_name': 'Decomposing wastes or polluting substances',
        'functional_role': 'producer',
        'confidence': 0.80,
        'pathway_context': 'Catechol meta-cleavage (xenobiotic degradation)',
    },

    # Photosynthesis → Climate regulation
    'M00161': {  # PSII
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.85,
        'pathway_context': 'Photosystem II (oxygenic photosynthesis)',
    },
    'M00163': {  # PSI
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.85,
        'pathway_context': 'Photosystem I (oxygenic photosynthesis)',
    },

    # Vitamin biosynthesis → Water quality (nutrient cycling)
    'M00125': {  # Riboflavin
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'producer',
        'confidence': 0.55,
        'pathway_context': 'Riboflavin biosynthesis',
    },
    'M00120': {  # B12
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'producer',
        'confidence': 0.60,
        'pathway_context': 'Cobalamin (B12) biosynthesis',
    },

    # Hydrogen metabolism → Climate regulation
    'M00627': {  # NiFe hydrogenase
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'transformer',
        'confidence': 0.60,
        'pathway_context': 'H2 uptake hydrogenase',
    },
    'M00628': {  # FeFe hydrogenase
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'transformer',
        'confidence': 0.60,
        'pathway_context': 'H2 production hydrogenase',
    },
}

# KEGG module definitions (step-based, from danaseq kegg_module_completeness.py)
KEGG_MODULES = {
    'M00175': ('Nitrogen fixation', [
        {'K02588'}, {'K02586'}, {'K02591'},
    ]),
    'M00528': ('Nitrification', [
        {'K10944', 'K10945', 'K10946'}, {'K10535'}, {'K00370', 'K00371'},
    ]),
    'M00529': ('Denitrification', [
        {'K00370', 'K00371', 'K02567', 'K02568'},
        {'K00368', 'K15864'},
        {'K04561', 'K02305'},
        {'K00376'},
    ]),
    'M00530': ('DNRA', [
        {'K00370', 'K00371', 'K02567', 'K02568'},
        {'K00362', 'K00363', 'K03385', 'K15876'},
    ]),
    'M00804': ('Anammox', [
        {'K20932', 'K20933', 'K20934'}, {'K20935'},
    ]),
    'M00531': ('Assimilatory nitrate reduction', [
        {'K00367', 'K10534'}, {'K00372', 'K00360'},
    ]),
    'M00595': ('Thiosulfate oxidation Sox', [
        {'K17222', 'K17223'}, {'K17224'}, {'K22622', 'K17225', 'K17226'},
    ]),
    'M00596': ('Dissimilatory sulfate reduction', [
        {'K00394', 'K00395'}, {'K11180', 'K11181'}, {'K00958'},
    ]),
    'M00176': ('Assimilatory sulfate reduction', [
        {'K00958', 'K00394'}, {'K00860'}, {'K00390'}, {'K00380', 'K00381'},
    ]),
    'M00567': ('Methanogenesis', [
        {'K00200', 'K00201', 'K00202', 'K00203'}, {'K00672'},
        {'K01499'}, {'K00319', 'K00320'},
        {'K00577', 'K00578', 'K00579', 'K00580', 'K00581', 'K00582', 'K00583', 'K00584'},
        {'K00399', 'K00401', 'K00402'},
    ]),
    'M00174': ('Methane oxidation', [
        {'K10944', 'K10945', 'K10946'},
        {'K16157', 'K16158', 'K16159', 'K16160', 'K16161', 'K16162'},
    ]),
    'M00346': ('Formaldehyde assimilation serine', [
        {'K00600'}, {'K00830'}, {'K00018'}, {'K01689'}, {'K01595'},
    ]),
    'M00165': ('Calvin cycle', [
        {'K00855'}, {'K01601', 'K01602'}, {'K00927'}, {'K00134', 'K00150'},
        {'K01623', 'K01624'}, {'K03841', 'K11532'}, {'K00615'},
        {'K01807', 'K01808'},
    ]),
    'M00173': ('rTCA', [
        {'K00169', 'K00170', 'K00171', 'K00172'}, {'K01648'},
        {'K15230', 'K15231'}, {'K00241', 'K00242'}, {'K01902', 'K01903'},
        {'K00174', 'K00175'},
    ]),
    'M00377': ('Wood-Ljungdahl', [
        {'K00198'}, {'K14138', 'K00197', 'K00194'}, {'K01938'},
        {'K01491'}, {'K00297'}, {'K15022', 'K15023'},
    ]),
    'M00376': ('3-HP', [
        {'K09709'}, {'K14468', 'K14469'}, {'K08691'},
        {'K01964', 'K01961', 'K01962'},
    ]),
    'M00548': ('Benzene degradation', [
        {'K16249', 'K16242', 'K16243', 'K16244'},
    ]),
    'M00568': ('Catechol ortho-cleavage', [
        {'K03381'}, {'K01856'}, {'K03464'},
    ]),
    'M00569': ('Catechol meta-cleavage', [
        {'K00446', 'K07104'}, {'K01666'},
    ]),
    'M00161': ('PSII', [
        {'K02703', 'K02706'}, {'K02705', 'K02704'},
    ]),
    'M00163': ('PSI', [
        {'K02689', 'K02690'}, {'K02694'},
    ]),
    'M00125': ('Riboflavin biosynthesis', [
        {'K01497'}, {'K14652'}, {'K00794'}, {'K00793'},
    ]),
    'M00120': ('Cobalamin B12 biosynthesis', [
        {'K02227'}, {'K02229', 'K02228'}, {'K02232', 'K02231'},
        {'K02225'}, {'K02233'},
    ]),
    'M00627': ('NiFe hydrogenase', [
        {'K06281', 'K06282'},
    ]),
    'M00628': ('FeFe hydrogenase', [
        {'K00533', 'K00534'},
    ]),
}

# ---------------------------------------------------------------------------
# FOAM L1 → CICES mapping
# Maps FOAM functional categories to ES codes for broader coverage.
# ---------------------------------------------------------------------------

FOAM_L1_ES = {
    '11_Nitrogen cycle': {
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'transformer',
        'confidence': 0.55,
    },
    '14_Methanogenesis': {
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'inhibitor',
        'confidence': 0.55,
    },
    '15_Methylotrophy': {
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'producer',
        'confidence': 0.50,
    },
    '18_Sulfur compounds metabolism': {
        'es_code': '2.3.5.1',
        'es_name': 'Controlling the chemical quality of freshwater',
        'functional_role': 'transformer',
        'confidence': 0.50,
    },
    '08_Hydrocarbon degradation': {
        'es_code': '2.1.1.1',
        'es_name': 'Decomposing wastes or polluting substances',
        'functional_role': 'producer',
        'confidence': 0.60,
    },
    '13_Hydrogen metabolism': {
        'es_code': '2.3.6.1',
        'es_name': 'Regulating our global climate',
        'functional_role': 'transformer',
        'confidence': 0.50,
    },
}


def load_foam(foam_path):
    """Load FOAM ontology: returns {KO: {L1, L2, L3, L4, Function, EC}}."""
    foam_kos = {}
    with open(foam_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ko = row.get('KO', '').strip()
            if ko and ko.startswith('K'):
                foam_kos[ko] = {
                    'L1': row.get('L1', '').strip(),
                    'L2': row.get('L2', '').strip(),
                    'L3': row.get('L3', '').strip(),
                    'L4': row.get('L4', '').strip(),
                    'Function': row.get('Function', '').strip(),
                    'EC': row.get('EC', '').strip(),
                }
    return foam_kos


def load_ontology_codes(ontology_path):
    """Load valid ES codes from ontology TSV."""
    codes = set()
    with open(ontology_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            codes.add(row['es_code'])
    return codes


def generate_kegg_rows():
    """Generate mapping rows from KEGG module definitions."""
    rows = []
    for module_id, es_info in KEGG_MODULE_ES.items():
        if module_id not in KEGG_MODULES:
            continue
        module_name, steps = KEGG_MODULES[module_id]
        all_kos = set()
        for step in steps:
            all_kos.update(step)

        for ko in sorted(all_kos):
            rows.append({
                'es_code': es_info['es_code'],
                'es_name': es_info['es_name'],
                'gene_id': ko,
                'gene_id_type': 'KO',
                'confidence': es_info['confidence'],
                'evidence_source': 'kegg_decoder',
                'functional_role': es_info['functional_role'],
                'pathway_context': es_info['pathway_context'],
                'notes': f'{module_name} ({module_id})',
            })
    return rows


def generate_foam_rows(foam_path):
    """Generate mapping rows from FOAM ontology."""
    foam_kos = load_foam(foam_path)
    rows = []
    for ko, info in foam_kos.items():
        l1 = info['L1']
        if l1 in FOAM_L1_ES:
            es_info = FOAM_L1_ES[l1]
            func_desc = info['Function'].split(';')[0].strip() if info['Function'] else ''
            rows.append({
                'es_code': es_info['es_code'],
                'es_name': es_info['es_name'],
                'gene_id': ko,
                'gene_id_type': 'KO',
                'confidence': es_info['confidence'],
                'evidence_source': 'foam',
                'functional_role': es_info['functional_role'],
                'pathway_context': f"{l1} / {info['L2']}" if info['L2'] else l1,
                'notes': func_desc,
            })
    return rows


def merge_and_deduplicate(kegg_rows, foam_rows):
    """Merge rows, boosting confidence where both sources agree."""
    # Index by (gene_id, es_code)
    merged = {}
    for row in kegg_rows:
        key = (row['gene_id'], row['es_code'])
        merged[key] = row.copy()

    for row in foam_rows:
        key = (row['gene_id'], row['es_code'])
        if key in merged:
            # Both sources agree — boost confidence
            existing = merged[key]
            boosted = min(0.95, max(existing['confidence'], row['confidence']) + 0.10)
            existing['confidence'] = boosted
            existing['evidence_source'] = 'kegg_decoder+foam'
            if row['notes'] and row['notes'] not in existing['notes']:
                existing['notes'] += f'; {row["notes"]}'
        else:
            merged[key] = row.copy()

    return sorted(merged.values(), key=lambda r: (r['es_code'], r['gene_id']))


def write_mapping(rows, output_path, valid_codes):
    """Write mapping TSV, filtering to valid ontology codes."""
    cols = ['es_code', 'es_name', 'gene_id', 'gene_id_type', 'confidence',
            'evidence_source', 'functional_role', 'pathway_context', 'notes']

    filtered = [r for r in rows if r['es_code'] in valid_codes]
    skipped = len(rows) - len(filtered)

    with open(output_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for row in filtered:
            f.write('\t'.join(str(row.get(c, '')) for c in cols) + '\n')

    return len(filtered), skipped


def main():
    parser = argparse.ArgumentParser(description='Bootstrap ES gene mapping table')
    parser.add_argument('--foam', required=True, help='FOAM ontology TSV')
    parser.add_argument('--ontology', required=True, help='ES ontology TSV')
    parser.add_argument('--output', required=True, help='Output mapping TSV')
    parser.add_argument('--kegg-only', action='store_true', help='Only KEGG, skip FOAM')
    args = parser.parse_args()

    valid_codes = load_ontology_codes(args.ontology)
    kegg_rows = generate_kegg_rows()
    print(f'KEGG-Decoder: {len(kegg_rows)} gene-ES mappings from {len(KEGG_MODULE_ES)} modules')

    if args.kegg_only:
        foam_rows = []
    else:
        foam_rows = generate_foam_rows(args.foam)
        print(f'FOAM: {len(foam_rows)} gene-ES mappings')

    merged = merge_and_deduplicate(kegg_rows, foam_rows)
    print(f'After merge/dedup: {len(merged)} mappings')

    written, skipped = write_mapping(merged, args.output, valid_codes)
    print(f'Written: {written} mappings (skipped {skipped} with invalid ES codes)')


if __name__ == '__main__':
    main()
