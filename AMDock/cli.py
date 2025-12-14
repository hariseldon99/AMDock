#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import subprocess
import sys

from AMDock.tools import PROJECT, BASE, Fix_PQR
from AMDock.variables import Variables


def load_defaults(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def copy_to_input(src, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    dst = os.path.join(dest_dir, os.path.basename(src))
    shutil.copy(src, dst)
    return dst


def run_cmd(cmd, cwd=None, capture=False):
    print('RUN:', ' '.join(cmd))
    res = subprocess.run(cmd, cwd=cwd, stdout=(subprocess.PIPE if capture else None), stderr=subprocess.STDOUT)
    return res.returncode


def prepare_inputs(v, proj, target_path, ligand_path, offtarget_path, args):
    # copy files into project input
    tgt = BASE()
    lig = BASE()
    tgt_path = copy_to_input(target_path, proj.input)
    tgt.get_data([tgt_path])
    tgt.input = tgt_path
    lig_path = copy_to_input(ligand_path, proj.input)
    lig.get_data([lig_path])
    lig.input = lig_path
    off = None
    if offtarget_path:
        off = BASE()
        off_path = copy_to_input(offtarget_path, proj.input)
        off.get_data([off_path])
        off.input = off_path

    # Prepare receptor: call pdb2pqr, Fix_PQR, prepare_receptor4
    if tgt.prepare:
        # pdb2pqr
        pdb2pqr = v.pdb2pqr_py
        pqr_out = os.path.join(proj.input, tgt.pqr)
        cmd = [pdb2pqr, '--titration-state-method=propka', '--noopt', '--drop-water', '--keep-chain',
               '--with-ph', str(args.get('pH', 7.4)), '--ff=%s' % args.get('forcefield', 'amber'), tgt.input, pqr_out]
        run_cmd(cmd, cwd=proj.input)
        # Fix_PQR callable
        metals = args.get('metals', [])
        Fix_PQR(tgt.input, pqr_out, metals)
        # prepare_receptor4
        prep_receptor = v.prepare_receptor4_py
        args_prep = [sys.executable, prep_receptor, '-r', tgt.pdb, '-v', '-U', 'nphs_lps_waters_nonstdres_deleteAltB']
        if args.get('metals_text'):
            args_prep += ['-p', args.get('metals_text')]
        run_cmd(args_prep, cwd=proj.input)

    # Prepare ligand
    if lig.prepare:
        prep_lig = v.prepare_ligand4_py
        lig_args = [sys.executable, prep_lig, '-l']
        if args.get('protonation') and args.get('protonation_program') == 'obabel':
            # run obabel protonation
            obabel = v.openbabel
            run_cmd([obabel, '-i', 'pdb', lig.input, '-opdb', '-O', lig.pdb, '-h', '-p', str(args.get('pH', 7.4))], cwd=proj.input)
            lig_args += [lig.pdb, '-v']
        else:
            lig_args += [lig.input, '-v']
        run_cmd(lig_args, cwd=proj.input)

    return tgt, lig, off


def prepare_grid_and_maps(v, proj, tgt, lig, off, args):
    # Prepare GPF and run AutoGrid (example for AutoDock4)
    prepare_gpf = v.prepare_gpf4_py
    spacing = args.get('spacing_autodock', v.spacing_autodock)
    # example prepare_gpf command
    gpf_args = [sys.executable, prepare_gpf, '-l', lig.pdbqt, '-r', tgt.pdbqt, '-p', 'spacing=%.3f' % spacing]
    run_cmd(gpf_args, cwd=proj.input)
    # run autogrid
    run_cmd([v.autogrid, '-p', tgt.auto_lig], cwd=proj.input)


def run_docking(v, proj, tgt, lig, off, args):
    prog = args.get('docking_program', 'AutoDock Vina')
    if 'Vina' in prog:
        # build vina command
        center = args.get('center')
        size = args.get('size')
        # example: vina --receptor receptor.pdbqt --ligand ligand.pdbqt --center_x ... --size_x ... --out out.pdbqt
        cmd = [v.vina_exec, '--receptor', tgt.pdbqt, '--ligand', lig.pdbqt]
        if center and size:
            cmd += ['--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
                    '--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2])]
        outp = os.path.join(proj.results, lig.name + '_vina_out.pdbqt')
        os.makedirs(proj.results, exist_ok=True)
        cmd += ['--out', outp]
        run_cmd(cmd, cwd=proj.input)
    else:
        # AutoDock4: prepare dpf and run autodock
        prepare_dpf = v.prepare_dpf_py
        dpf_args = [sys.executable, prepare_dpf, '-l', lig.pdbqt, '-r', tgt.pdbqt, '-p', 'spacing=%.3f' % v.spacing_autodock]
        run_cmd(dpf_args, cwd=proj.input)
        run_cmd([v.autodock, '-p', tgt.dpf], cwd=proj.input)


def main():
    parser = argparse.ArgumentParser(description='AMDock command-line runner (partial feature set)')
    parser.add_argument('--config', help='JSON defaults file', default=os.path.join(os.path.dirname(__file__), 'data', 'cli_defaults.json'))
    parser.add_argument('--wdir', help='Project location', default='.')
    parser.add_argument('--project', help='Project name', default='Project_Docking')
    parser.add_argument('--target', help='Target protein file', required=True)
    parser.add_argument('--ligand', help='Ligand file', required=True)
    parser.add_argument('--offtarget', help='Off-target protein file', default=None)
    parser.add_argument('--steps', help='Comma-separated steps: prepare,grid,dock', default='prepare,grid,dock')
    parser.add_argument('--docking_program', help='Docking program (AutoDock Vina or AutoDock4)', default=None)
    args = parser.parse_args()

    defaults = load_defaults(args.config)
    v = Variables()
    # build project
    proj = PROJECT()
    proj.get_loc([args.wdir])
    proj.get_info(args.project)
    os.makedirs(proj.input, exist_ok=True)
    os.makedirs(proj.results, exist_ok=True)

    # merge CLI args with defaults
    merged = dict(defaults)
    merged.update({k: v for k, v in vars(args).items() if v is not None})

    tgt, lig, off = prepare_inputs(v, proj, args.target, args.ligand, args.offtarget, merged)
    steps = [s.strip() for s in args.steps.split(',')]
    if 'grid' in steps:
        prepare_grid_and_maps(v, proj, tgt, lig, off, merged)
    if 'dock' in steps:
        run_docking(v, proj, tgt, lig, off, merged)


if __name__ == '__main__':
    main()
