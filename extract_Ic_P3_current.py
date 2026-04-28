# -*- coding: utf-8 -*-
from __future__ import print_function

"""
Extract Ic on the manually inserted body cohesive path from an Abaqus ODB.

This version is adapted for the current Chapter 4 models:
- interface uses surface-based cohesive contact, not cohesive elements
- the body path cohesive seam is inserted manually in Mesh Edit
- many times there is NO usable element set for the inserted seam in the ODB

So the script uses the following priority:
1. try a matching element set if one exists
2. otherwise collect cohesive elements and filter them by x-coordinate of the path
3. if the x-filter fails, fall back to all cohesive elements and warn clearly

Default target = P3, x = 0.25 mm

Usage examples
- abaqus python extract_Ic_P3_current.py
- abaqus python extract_Ic_P3_current.py Job-1.odb
- abaqus python extract_Ic_P3_current.py Job-1.odb cooling
- abaqus python extract_Ic_P3_current.py Job-1.odb cooling 0.25
- abaqus python extract_Ic_P3_current.py Job-1.odb cooling 0.25 0.002

Outputs
- Ic_P3_history.csv
- Ic_P3_summary.txt
"""

import os
import sys
import csv

try:
    from odbAccess import openOdb
except Exception:
    raise RuntimeError('This script must be run with Abaqus Python, e.g. "abaqus python extract_Ic_P3_current.py"')

DEFAULT_ODB = 'Job-1.odb'
DEFAULT_STEP = 'cooling'
DEFAULT_X_TARGET = 0.25
DEFAULT_X_TOL = 0.002
TARGET_SET_CANDIDATES = ('SET_VERTICAL_PATH_P3', 'SET_VERTICAL_PATH', 'VERTICAL_PATH_P3', 'VERTICAL_PATH')
COH_TYPES = ('COH2D4', 'COH2D4P', 'COHAX4', 'COH3D8', 'COH3D6', 'COH3D4')


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
    x_target = DEFAULT_X_TARGET
    x_tol = DEFAULT_X_TOL

    non_odb = []
    for a in clean:
        if a.lower().endswith('.odb'):
            odb_path = a
        else:
            non_odb.append(a)

    if len(non_odb) >= 1:
        step_name = non_odb[0]
    if len(non_odb) >= 2:
        try:
            x_target = float(non_odb[1])
        except Exception:
            pass
    if len(non_odb) >= 3:
        try:
            x_tol = float(non_odb[2])
        except Exception:
            pass

    return odb_path, step_name, x_target, x_tol, clean


def find_step_name(odb, preferred_name):
    if preferred_name in odb.steps.keys():
        return preferred_name
    upper_map = _upper_keys(odb.steps)
    if preferred_name.upper() in upper_map:
        return upper_map[preferred_name.upper()]
    return list(odb.steps.keys())[-1]


def find_region_set(odb, candidate_names):
    for target_name in candidate_names:
        target_u = target_name.upper()

        for inst_name, inst in odb.rootAssembly.instances.items():
            upper_map = _upper_keys(inst.elementSets)
            if target_u in upper_map:
                real_name = upper_map[target_u]
                return inst.elementSets[real_name], 'instance elementSet %s:%s' % (inst_name, real_name)

        upper_map = _upper_keys(odb.rootAssembly.elementSets)
        if target_u in upper_map:
            real_name = upper_map[target_u]
            return odb.rootAssembly.elementSets[real_name], 'assembly elementSet %s' % real_name

    return None, None


def collect_all_cohesive_elements(odb):
    inst_to_elems = {}
    total = 0
    for inst_name, inst in odb.rootAssembly.instances.items():
        elems = []
        for e in inst.elements:
            etype = str(e.type).upper()
            if etype in COH_TYPES or etype.startswith('COH'):
                elems.append(e)
        if elems:
            inst_to_elems[inst_name] = elems
            total += len(elems)
    return inst_to_elems, total


def build_node_coord_map(inst):
    d = {}
    for n in inst.nodes:
        d[n.label] = n.coordinates
    return d


def element_centroid_x(elem, node_coord_map):
    xs = []
    for nlab in elem.connectivity:
        if nlab in node_coord_map:
            xs.append(node_coord_map[nlab][0])
    if not xs:
        return None
    return sum(xs) / float(len(xs))


def collect_cohesive_labels_near_x(odb, x_target, x_tol):
    label_map = {}
    total = 0
    bbox = None

    for inst_name, inst in odb.rootAssembly.instances.items():
        node_coord_map = build_node_coord_map(inst)
        labels = set()
        xs_all = []
        ys_all = []

        for e in inst.elements:
            etype = str(e.type).upper()
            if not (etype in COH_TYPES or etype.startswith('COH')):
                continue

            xc = element_centroid_x(e, node_coord_map)
            if xc is None:
                continue
            if abs(xc - x_target) <= x_tol:
                labels.add(e.label)
                xs_all.append(xc)
                for nlab in e.connectivity:
                    if nlab in node_coord_map:
                        ys_all.append(node_coord_map[nlab][1])

        if labels:
            label_map[inst_name] = labels
            total += len(labels)
            if xs_all and ys_all:
                x_min = min(xs_all)
                x_max = max(xs_all)
                y_min = min(ys_all)
                y_max = max(ys_all)
                bbox = (x_min, x_max, y_min, y_max, inst_name)

    return label_map, total, bbox


def get_scalar_field(frame):
    if 'MAXSCRT' in frame.fieldOutputs.keys():
        return frame.fieldOutputs['MAXSCRT'], 'MAXSCRT'

    if 'DMICRT' in frame.fieldOutputs.keys():
        fo = frame.fieldOutputs['DMICRT']
        try:
            sfield = fo.getScalarField(componentLabel='MAXSCRT')
            return sfield, 'DMICRT.MAXSCRT'
        except Exception:
            return fo, 'DMICRT'

    raise RuntimeError('Neither MAXSCRT nor DMICRT is available in the ODB output.')


def scalarize_data(vdata):
    if isinstance(vdata, (float, int)):
        return float(vdata)
    try:
        if len(vdata) == 1:
            return float(vdata[0])
        return float(vdata[0])
    except Exception:
        return float(vdata)


def write_csv(csv_path, rows):
    with open(csv_path, 'w') as f:
        w = csv.writer(f)
        w.writerow(['frame_index', 'frame_value', 'Ic_frame_max', 'instance_name', 'element_label'])
        for r in rows:
            w.writerow(r)


def write_summary(txt_path, odb_path, step_name, source_name, method_desc, x_target, x_tol, bbox, rows):
    valid = [r for r in rows if r[2] is not None]
    with open(txt_path, 'w') as f:
        f.write('ODB: %s\n' % odb_path)
        f.write('Step: %s\n' % step_name)
        f.write('Field source: %s\n' % source_name)
        f.write('Selection method: %s\n' % method_desc)
        f.write('Target x: %.12g\n' % x_target)
        f.write('x tolerance: %.12g\n' % x_tol)
        if bbox is not None:
            f.write('Selected cohesive path bbox: x=[%.12g, %.12g], y=[%.12g, %.12g], instance=%s\n'
                    % (bbox[0], bbox[1], bbox[2], bbox[3], bbox[4]))
        f.write('Number of frames: %d\n' % len(rows))

        if not valid:
            f.write('No valid Ic data found.\n')
            return

        max_row = max(valid, key=lambda x: x[2])
        f.write('Ic_overall_max: %.12g\n' % max_row[2])
        f.write('At frame_index: %d\n' % max_row[0])
        f.write('At frame_value: %.12g\n' % max_row[1])
        f.write('At instance: %s\n' % max_row[3])
        f.write('At element_label: %s\n' % str(max_row[4]))
        if max_row[2] >= 1.0:
            f.write('Interpretation: Ic has reached initiation threshold on this path.\n')
        else:
            f.write('Interpretation: Ic has not yet reached initiation threshold on this path.\n')


def main():
    odb_path, step_name, x_target, x_tol, clean_args = parse_user_args(sys.argv)

    if not os.path.exists(odb_path):
        raise RuntimeError(
            'ODB file not found: %s\n'
            'Clean args seen by script: %s\n'
            'Use for example: abaqus python extract_Ic_P3_current.py Job-1.odb cooling 0.25 0.002'
            % (odb_path, clean_args)
        )

    odb = openOdb(path=odb_path, readOnly=True)

    try:
        step_name = find_step_name(odb, step_name)
        step = odb.steps[step_name]

        region_set, region_desc = find_region_set(odb, TARGET_SET_CANDIDATES)
        use_region_set = False
        use_label_filter = False
        label_map = {}
        method_desc = ''
        bbox = None

        if region_set is not None:
            use_region_set = True
            method_desc = region_desc
            print('Using region: %s' % region_desc)
        else:
            label_map, total, bbox = collect_cohesive_labels_near_x(odb, x_target, x_tol)
            if total > 0:
                use_label_filter = True
                method_desc = 'cohesive elements filtered by centroid x near %.12g +/- %.12g' % (x_target, x_tol)
                print('Using x-filtered cohesive elements: count = %d' % total)
                if bbox is not None:
                    print('Selected path bbox: x=[%.6f, %.6f], y=[%.6f, %.6f], instance=%s'
                          % (bbox[0], bbox[1], bbox[2], bbox[3], bbox[4]))
            else:
                all_coh, total_all = collect_all_cohesive_elements(odb)
                if total_all == 0:
                    raise RuntimeError('No cohesive elements were found in the ODB. Check the manual seam insertion.')
                label_map = {}
                for inst_name, elems in all_coh.items():
                    label_map[inst_name] = set([e.label for e in elems])
                use_label_filter = True
                method_desc = 'WARNING: fallback to ALL cohesive elements in the ODB'
                print('Warning: no path set and no x-filtered cohesive elements were found.')
                print('Fallback to ALL cohesive elements in the ODB. Verify that the model contains only one body cohesive seam.')

        rows = []
        source_name_used = None

        for i, frame in enumerate(step.frames):
            scalar_field, source_name = get_scalar_field(frame)
            source_name_used = source_name

            if use_region_set:
                subset = scalar_field.getSubset(region=region_set)
                values = subset.values
            else:
                values = []
                for v in scalar_field.values:
                    try:
                        inst_name = v.instance.name
                    except Exception:
                        continue
                    if inst_name in label_map and v.elementLabel in label_map[inst_name]:
                        values.append(v)

            if not values:
                rows.append([i, frame.frameValue, None, '', ''])
                continue

            best_val = None
            best_inst = ''
            best_label = ''

            for v in values:
                sval = scalarize_data(v.data)
                if best_val is None or sval > best_val:
                    best_val = sval
                    try:
                        best_inst = v.instance.name
                    except Exception:
                        best_inst = ''
                    best_label = getattr(v, 'elementLabel', '')

            rows.append([i, frame.frameValue, best_val, best_inst, best_label])

        csv_path = 'Ic_P3_history.csv'
        txt_path = 'Ic_P3_summary.txt'
        write_csv(csv_path, rows)
        write_summary(txt_path, odb_path, step_name, source_name_used, method_desc, x_target, x_tol, bbox, rows)

        valid = [r for r in rows if r[2] is not None]
        if valid:
            max_row = max(valid, key=lambda x: x[2])
            print('Ic extraction finished.')
            print('Step used: %s' % step_name)
            print('Field source: %s' % source_name_used)
            print('Selection method: %s' % method_desc)
            print('Ic_overall_max = %.12g' % max_row[2])
            print('frame_index = %d, frame_value = %.12g' % (max_row[0], max_row[1]))
            if max_row[2] >= 1.0:
                print('Interpretation: Ic >= 1.0, the body path has reached initiation.')
            else:
                print('Interpretation: Ic < 1.0, the body path is activated but has not yet initiated.')
            print('CSV written to: %s' % csv_path)
            print('Summary written to: %s' % txt_path)
        else:
            print('No valid Ic values were found. Check MAXSCRT field output and the manual seam insertion.')

    finally:
        odb.close()


if __name__ == '__main__':
    main()
