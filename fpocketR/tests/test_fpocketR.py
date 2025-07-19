"""
Unit and regression test for the fpocketR package.
"""

import sys
import subprocess
from pathlib import Path
import os
import pytest
import pandas as pd
import atexit


# --- Cleanup orphaned fpocket output files after tests ---
@pytest.fixture(scope="session", autouse=True)
def cleanup_fpocket_files():
    files_to_remove = [
        Path(__file__).parent.parent.parent / "2l1v_clean.pdb",
        Path(__file__).parent.parent.parent / "8f4o_clean.pdb"
    ]
    def remove_files():
        for f in files_to_remove:
            if f.exists():
                f.unlink()
    atexit.register(remove_files)


# --- Basic Import Test ---
def test_fpocketR_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "fpocketR" in sys.modules



# --- Helper Functions ---
def get_multistate_subdirs(multistate_dir):
    """Return sorted list of 2l1v_clean_stateX_out subdirectories."""
    return sorted([d for d in multistate_dir.iterdir() if d.is_dir() and d.name.startswith("2l1v_clean_state")])

def get_expected_state_dirnames():
    return [f"2l1v_clean_state{i}_out" for i in range(1, 4)]

# --- CSV Comparison Helper ---
def tolerant_csv_compare(file1, file2, atol=10):
    """Compare two CSVs using pandas, allowing small differences in float columns."""
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    if df1.shape != df2.shape or list(df1.columns) != list(df2.columns):
        print(f"[DEBUG] CSV shape/columns mismatch: {file1} vs {file2}")
        return False
    # If 'Volume' column exists, allow higher tolerance for it
    if 'Volume' in df1.columns:
        try:
            # Compare all columns except 'Volume'
            cols = [c for c in df1.columns if c != 'Volume']
            pd.testing.assert_frame_equal(df1[cols], df2[cols], atol=atol, rtol=0, check_dtype=False, check_exact=False)
            # Compare 'Volume' with even higher tolerance
            pd.testing.assert_series_equal(df1['Volume'], df2['Volume'], atol=atol * 3, rtol=0, check_dtype=False, check_exact=False)
            return True
        except AssertionError as e:
            print(f"[DEBUG] CSV difference in {file1} vs {file2}:")
            print(e)
            return False
    else:
        try:
            pd.testing.assert_frame_equal(df1, df2, atol=atol, rtol=0, check_dtype=False, check_exact=False)
            return True
        except AssertionError as e:
            print(f"[DEBUG] CSV difference in {file1} vs {file2}:")
            print(e)
            return False


# --- Single-State Output Fixture (Alternate Command) ---
@pytest.fixture(scope="session")
def single_state_alt_output(tmp_path_factory):
    """Run fpocketR with alternate command and return output dir."""
    repo_root = Path(__file__).parent.parent.parent.resolve()
    data_dir = repo_root / "fpocketR" / "data"
    pdb_file = data_dir / "8f4o.pdb"
    al_file = data_dir / "2gdi.pdb"
    tmp_path = tmp_path_factory.mktemp("single_state_alt")
    out_dir = tmp_path / "8f4o_clean_out"
    out_dir.mkdir(parents=True, exist_ok=True)
    rel_out_dir = os.path.relpath(out_dir, repo_root)
    cmd = [
        "python", "-m", "fpocketR",
        "-pdb", str(pdb_file),
        "-al", str(al_file),
        "-l", "no",
        "-dpi", "50",
        "-o", rel_out_dir
    ]
    result = subprocess.run(cmd, check=True, cwd=str(repo_root), capture_output=True, text=True)
    return out_dir


# --- Alternate Single-State Output Tests ---
def test_single_state_alt_pdb_match_reference(single_state_alt_output):
    """Compare generated 8f4o_out_real_sphere.pdb to reference file."""
    out_dir = single_state_alt_output
    repo_root = Path(__file__).parent.parent.parent.resolve()
    data_dir = repo_root / "fpocketR" / "data"
    generated_pdbs = list(out_dir.rglob("8f4o_out_real_sphere.pdb"))
    assert generated_pdbs, f"No generated 8f4o_out_real_sphere.pdb file found in {out_dir}"
    generated_pdb = generated_pdbs[0]
    reference_pdb = data_dir / "8f4o_clean_out" / "8f4o_out_real_sphere.pdb"
    assert reference_pdb.exists(), f"Reference file not found: {reference_pdb}"
    with open(generated_pdb, "r") as f1, open(reference_pdb, "r") as f2:
        gen_lines = [line.strip() for line in f1 if line.strip()]
        ref_lines = [line.strip() for line in f2 if line.strip()]
    assert gen_lines == ref_lines, f"Generated PDB does not match reference PDB in {out_dir}"


def test_single_state_alt_3d_files_exist(single_state_alt_output):
    """Test that 3D files are generated for alternate single-state output."""
    out_dir = single_state_alt_output
    files = list(out_dir.rglob("*"))
    files = [f for f in files if f.is_file()]
    assert any(f.name.endswith("_3D_50.png") for f in files), f"No _3D_50.png in {out_dir}"
    assert any(f.name.endswith("_real_sphere.pdb") for f in files), f"No _real_sphere.pdb in {out_dir}"
    assert any(f.name.endswith("_real_sphere.pse") for f in files), f"No _real_sphere.pse in {out_dir}"


def test_single_state_alt_csv_match_reference(single_state_alt_output):
    """Compare generated *_pocket_characteristics.csv to reference file for alt single-state output."""
    out_dir = single_state_alt_output
    repo_root = Path(__file__).parent.parent.parent.resolve()
    data_dir = repo_root / "fpocketR" / "data"
    generated_csvs = list(out_dir.rglob("8f4o_out_pocket_characteristics.csv"))
    assert generated_csvs, f"No generated 8f4o_out_pocket_characteristics.csv file found in {out_dir}"
    generated_csv = generated_csvs[0]
    reference_csv = data_dir / "8f4o_clean_out" / "8f4o_out_pocket_characteristics.csv"
    assert reference_csv.exists(), f"Reference file not found: {reference_csv}"
    assert tolerant_csv_compare(generated_csv, reference_csv), f"Generated CSV does not match reference CSV in {out_dir}"


# --- Multistate Output Fixture ---

import pytest

@pytest.fixture(scope="session")
def multistate_output(tmp_path_factory):
    """Run the multistate fpocketR command ONCE and return output dir. Prints diagnostics and skips tests if fails."""
    tmp_path = tmp_path_factory.mktemp("multistate")
    repo_root = Path(__file__).parent.parent.parent.resolve()
    data_dir = repo_root / "fpocketR" / "data"
    pdb_file = (data_dir / "2l1v.pdb").resolve()
    ss_file = (data_dir / "2l1v.nsd").resolve()
    out_dir = tmp_path / "2l1v_multistate"
    out_dir.mkdir(parents=True, exist_ok=True)
    rel_out_dir = os.path.relpath(out_dir, repo_root)
    cmd = [
        "python", "-m", "fpocketR",
        "-pdb", str(pdb_file),
        "-ss", str(ss_file),
        "-s", "0", "-dpi", "50", "-o", rel_out_dir
    ]
    try:
        result = subprocess.run(cmd, check=True, cwd=str(repo_root), capture_output=True, text=True)
        return out_dir
    except subprocess.CalledProcessError as e:
        import sys
        sys.stderr.write("\n[ERROR] Multistate fpocketR subprocess failed!\n")
        sys.stderr.write(f"Command: {e.cmd}\n")
        sys.stderr.write(f"Return code: {e.returncode}\n")
        sys.stderr.write(f"Stdout:\n{e.output}\n")
        sys.stderr.write(f"Stderr:\n{e.stderr}\n")
        pytest.fail(f"Multistate fpocketR run failed: {e}")


# --- Multistate Output Tests ---
MULTISTATE_DIR = Path(__file__).parent.parent.parent / "fpocketR" / "data" / "2l1v_multistate"

@pytest.mark.parametrize("state_idx", range(1, 4))
def test_multistate_3d_files_exist(multistate_output, state_idx):
    """Test that 3D files exist in each state output directory."""
    out_dir = multistate_output
    gen_dir = out_dir / f"2l1v_clean_state{state_idx}_out"
    assert gen_dir.exists(), f"Missing generated subdir: {gen_dir}"
    files = [f for f in gen_dir.iterdir() if f.is_file()]
    assert any(f.name.endswith("_3D_50.png") for f in files), f"No _3D_50.png in {gen_dir}"
    assert any(f.name.endswith("_real_sphere.pdb") for f in files), f"No _real_sphere.pdb in {gen_dir}"
    assert any(f.name.endswith("_real_sphere.pse") for f in files), f"No _real_sphere.pse in {gen_dir}"


@pytest.mark.parametrize("state_idx", range(1, 4))
def test_multistate_2d_files_exist(multistate_output, state_idx):
    out_dir = multistate_output
    gen_dir = out_dir / f"2l1v_clean_state{state_idx}_out"
    assert gen_dir.exists(), f"Missing generated subdir: {gen_dir}"
    files = [f for f in gen_dir.iterdir() if f.is_file()]
    assert any(f.name.endswith("_2D.png") for f in files), f"No _2D.png in {gen_dir}"
    assert any(f.name.endswith("_2D.svg") for f in files), f"No _2D.svg in {gen_dir}"


@pytest.mark.parametrize("state_idx", range(1, 4))
def test_multistate_1d_csv_match_reference(multistate_output, state_idx):
    out_dir = multistate_output
    gen_dir = out_dir / f"2l1v_clean_state{state_idx}_out"
    ref_dir = MULTISTATE_DIR / f"2l1v_clean_state{state_idx}_out"
    gen_csvs = list(gen_dir.glob("*_pocket_characteristics.csv"))
    ref_csvs = list(ref_dir.glob("*_pocket_characteristics.csv"))
    assert gen_csvs, f"No generated *_pocket_characteristics.csv in {gen_dir}"
    assert ref_csvs, f"No reference *_pocket_characteristics.csv in {ref_dir}"
    gen_csv = gen_csvs[0]
    ref_csv = ref_csvs[0]
    assert tolerant_csv_compare(gen_csv, ref_csv), f"CSV mismatch in {gen_dir}"


def test_multistate_summary_files_exist(multistate_output):
    """Test that summary files exist in the output directory."""
    out_dir = multistate_output
    gen_files = [f for f in out_dir.iterdir() if f.is_file()]
    # Accept any file matching *_all_states_out_real_sphere.pse
    assert any(f.name.endswith("_all_states_out_real_sphere.pse") for f in gen_files), "Missing *_all_states_out_real_sphere.pse"
    assert any(f.name.endswith("_pocket_density.png") for f in gen_files), "Missing _pocket_density.png"
    assert any(f.name.endswith("all_states_pocket_characteristics.csv") for f in gen_files), "Missing all_states_pocket_characteristics.csv"