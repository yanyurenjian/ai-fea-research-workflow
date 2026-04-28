# -*- coding: utf-8 -*-
from __future__ import print_function

"""
Extract Id history from surface-based cohesive contact in an Abaqus ODB.

Id(frame) = max(interface contact initiation scalar on the whole interface)

Robustness improvements versus earlier version:
1. does not crash when the first few frames do not contain the contact field output
2. searches several candidate field names and also scans available keys
3. if the field is missing in some frames, records Id=0.0 for those frames and continues
4. keeps the same interface-surface selection logic used in the current Chapter 4 models

Usage
- abaqus python extract_Id_current_fix.py
- abaqus python extract_Id_current_fix.py Job-1.odb
- abaqus python extract_Id_current_fix.py Job-1.odb cooling

Outputs
- Id_history.csv
- Id_summary.txt
"""

import os
import sys
import csv

try:
    from odbAccess import openOdb
except Exception:
    raise RuntimeError('This script must be run with Abaqus Python, e.g. "abaqus python extract_Id_current_fix.py"')

DEFAULT_ODB = 'Job-1.odb'
DEFAULT_STEP = 'cooling'
SURFACE_CANDIDATES = (
    'SURF_SUB_WAVETOP',
    'SURF_COAT_BOTTOM',
    'SURF_SUB_WAVETOP_ALL',
    'SURF_COAT_BOTTOM_ALL',
)
FIELD_CANDIDATES = (
    'CSQUADSCRT',
    'CQUADSCRT',
)
MISSING_FIELD_VALUE = 0.0


def _upper_keys(d):
    return dict((k.upper(), k) for k in d.keys())


def parse_user_args(argv):
    clean = []
    for a in argv[1:]:
        if not a:
            continue
        if a.startswith('-'):
            continue
        clean.append(a)

    odb_path = DEFAULT_ODB
    step_name = DEFAULT_STEP

    non_odb = []
    for a in clean:
        if a.lower().endswith('.odb'):
            odb_path = a
        else:
            non_odb.append(a)

    if len(non_odb) >= 1:
        step_name = non_odb[0]

    return odb_path, step_name, clean


def find_step_name(odb, preferred_name):
    if preferred_name in odb.steps.keys():
        return preferred_name
    upper_map = _upper_keys(odb.steps)
    if preferred_name.upper() in upper_map:
        return upper_map[preferred_name.upper()]
    return list(odb.steps.keys())[-1]


def find_surface_region(odb, candidate_names):
    asm = odb.rootAssembly

    if hasattr(asm, 'surfaces'):
        asm_map = _upper_keys(asm.surfaces)
        for target_name in candidate_names:
            tu = target_name.upper()
            if tu in asm_map:
                rn = asm_map[tu]
                return asm.surfaces[rn], 'assembly surface %s' % rn

    for inst_name, inst in asm.instances.items():
        if hasattr(inst, 'surfaces'):
            inst_map = _upper_keys(inst.surfaces)
            for target_name in candidate_names:
                tu = target_name.upper()
                if tu in inst_map:
                    rn = inst_map[tu]
                    return inst.surfaces[rn], 'instance surface %s:%s' % (inst_name, rn)

    return None, None


def discover_field_name(frame):
    keys = frame.fieldOutputs.keys()
    upper_map = _upper_keys(frame.fieldOutputs)

    for name in FIELD_CANDIDATES:
        if name.upper() in upper_map:
            real_name = upper_map[name.upper()]
            return frame.fieldOutputs[real_name], real_name

    # fuzzy fallback: anything containing QUADSCRT and starting with C
    for k in keys:
        ku = k.upper()
        if ku.startswith('C') and 'QUADSCRT' in ku:
            return frame.fieldOutputs[k], k

    return None, None


def scalarize_data(vdata):
    if isinstance(vdata, (float, int)):
        return float(vdata)
    try:
        if len(vdata) == 1:
            return float(vdata[0])
        return float(vdata[0])
    except Exception:
        return float(vdata)


def best_value_from_values(values):
    best_val = None
    best_inst = ''
    best_node = ''
    best_elem = ''

    for v in values:
        try:
            sval = scalarize_data(v.data)
        except Exception:
            continue
        if best_val is None or sval > best_val:
            best_val = sval
            try:
                best_inst = v.instance.name
            except Exception:
                best_inst = ''
            best_node = getattr(v, 'nodeLabel', '')
            best_elem = getattr(v, 'elementLabel', '')

    return best_val, best_inst, best_node, best_elem


def write_csv(csv_path, rows):
    with open(csv_path, 'w') as f:
        w = csv.writer(f)
        w.writerow(['frame_index', 'frame_value', 'Id_frame_max', 'instance_name', 'node_label', 'element_label', 'field_name', 'status'])
        for r in rows:
            w.writerow(r)


def write_summary(txt_path, odb_path, step_name, field_name, method_desc, rows, missing_count, discovered_frame, available_keys_note):
    valid = [r for r in rows if r[2] is not None]
    with open(txt_path, 'w') as f:
        f.write('ODB: %s\n' % odb_path)
        f.write('Step: %s\n' % step_name)
        f.write('Field source used: %s\n' % (field_name if field_name else 'NOT FOUND'))
        f.write('Selection method: %s\n' % method_desc)
        f.write('Number of frames: %d\n' % len(rows))
        f.write('Frames without field output: %d\n' % missing_count)
        if discovered_frame is not None:
            f.write('First frame where field was found: %d\n' % discovered_frame)
        if available_keys_note:
            f.write('Field keys seen in last inspected frame: %s\n' % available_keys_note)

        if not valid:
            f.write('No valid Id data found.\n')
            return

        max_row = max(valid, key=lambda x: x[2])
        f.write('Id_overall_max: %.12g\n' % max_row[2])
        f.write('At frame_index: %d\n' % max_row[0])
        f.write('At frame_value: %.12g\n' % max_row[1])
        f.write('At instance: %s\n' % max_row[3])
        f.write('At node_label: %s\n' % str(max_row[4]))
        f.write('At element_label: %s\n' % str(max_row[5]))

        first_ge1 = None
        for r in valid:
            if r[2] >= 1.0:
                first_ge1 = r
                break
        if first_ge1 is None:
            f.write('Interpretation: Id has not yet reached initiation threshold on the interface.\n')
        else:
            f.write('First frame with Id >= 1: frame_index=%d, frame_value=%.12g\n' % (first_ge1[0], first_ge1[1]))
            f.write('Interpretation: Id has reached initiation threshold on the interface.\n')


def main():
    odb_path, step_name, clean_args = parse_user_args(sys.argv)

    if not os.path.exists(odb_path):
        raise RuntimeError(
            'ODB file not found: %s\n'
            'Clean args seen by script: %s\n'
            'Use for example: abaqus python extract_Id_current_fix.py Job-1.odb cooling'
            % (odb_path, clean_args)
        )

    odb = openOdb(path=odb_path, readOnly=True)

    try:
        step_name = find_step_name(odb, step_name)
        step = odb.steps[step_name]

        surf_region, region_desc = find_surface_region(odb, SURFACE_CANDIDATES)
        use_region = surf_region is not None
        if use_region:
            method_desc = region_desc
            print('Using interface region: %s' % region_desc)
        else:
            method_desc = 'WARNING: fallback to all available interface initiation field values in each frame'
            print('Warning: named interface surface was not found in ODB.')
            print('Fallback to all available interface initiation field values in each frame.')

        rows = []
        field_name_used = None
        missing_count = 0
        discovered_frame = None
        available_keys_note = ''

        for i, frame in enumerate(step.frames):
            scalar_field, field_name = discover_field_name(frame)
            available_keys_note = ', '.join(sorted(frame.fieldOutputs.keys()))

            if scalar_field is None:
                missing_count += 1
                rows.append([i, frame.frameValue, MISSING_FIELD_VALUE, '', '', '', '', 'field_missing_assumed_zero'])
                continue

            if field_name_used is None:
                field_name_used = field_name
                discovered_frame = i
                print('Detected interface initiation field: %s (first found at frame %d)' % (field_name_used, i))

            if use_region:
                try:
                    subset = scalar_field.getSubset(region=surf_region)
                    values = subset.values
                except Exception:
                    values = scalar_field.values
                    if 'getsubset failed' not in method_desc.lower():
                        method_desc = method_desc + ' ; getSubset failed, fallback to all available field values'
            else:
                values = scalar_field.values

            if not values:
                rows.append([i, frame.frameValue, MISSING_FIELD_VALUE, '', '', '', field_name, 'empty_values_assumed_zero'])
                continue

            best_val, best_inst, best_node, best_elem = best_value_from_values(values)
            if best_val is None:
                best_val = MISSING_FIELD_VALUE
                rows.append([i, frame.frameValue, best_val, '', '', '', field_name, 'unreadable_values_assumed_zero'])
            else:
                rows.append([i, frame.frameValue, best_val, best_inst, best_node, best_elem, field_name, 'ok'])

        csv_path = 'Id_history.csv'
        txt_path = 'Id_summary.txt'
        write_csv(csv_path, rows)
        write_summary(txt_path, odb_path, step_name, field_name_used, method_desc, rows, missing_count, discovered_frame, available_keys_note)

        valid = [r for r in rows if r[2] is not None]
        if valid:
            max_row = max(valid, key=lambda x: x[2])
            first_ge1 = None
            for r in valid:
                if r[2] >= 1.0:
                    first_ge1 = r
                    break
            print('Id extraction finished.')
            print('Step used: %s' % step_name)
            print('Field source used: %s' % (field_name_used if field_name_used else 'NOT FOUND'))
            print('Selection method: %s' % method_desc)
            print('Frames without field output = %d' % missing_count)
            print('Id_overall_max = %.12g' % max_row[2])
            print('frame_index = %d, frame_value = %.12g' % (max_row[0], max_row[1]))
            if first_ge1 is None:
                print('Interpretation: Id < 1.0 for all frames, the interface has not yet initiated.')
            else:
                print('First frame with Id >= 1.0: frame_index = %d, frame_value = %.12g' % (first_ge1[0], first_ge1[1]))
                print('Interpretation: Id has reached initiation on the interface.')
            print('CSV written to: %s' % csv_path)
            print('Summary written to: %s' % txt_path)
        else:
            print('No valid Id values were found. Check interface field output requests and cohesive contact definition.')

    finally:
        odb.close()


if __name__ == '__main__':
    main()
