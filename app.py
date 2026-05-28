from flask import Flask, request, jsonify, render_template
import os
import subprocess
import base64
import uuid
import sys
import csv
from datetime import datetime

# ─────────────────────────────────────────────────────────────────
# app.py — the Flask web server
#
# Think of this file as the "manager" of your restaurant.
# It receives files from the browser, tells the kitchen (your scripts)
# to run, then sends results back to the browser.
# ─────────────────────────────────────────────────────────────────

app = Flask(__name__)

# Where this file lives (project root)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Folders for temporary uploads and results
UPLOAD_FOLDER  = os.path.join(BASE_DIR, "data", "uploads")
RESULTS_FOLDER = os.path.join(BASE_DIR, "data", "results")
IMAGES_FOLDER  = os.path.join(BASE_DIR, "data", "images")
SCRIPTS_FOLDER = os.path.join(BASE_DIR, "scripts")

# Create these folders if they don't exist yet
os.makedirs(UPLOAD_FOLDER,  exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)
os.makedirs(IMAGES_FOLDER,  exist_ok=True)


# ── Route 1: Serve the main webpage ───────────────────────────────
# When someone visits your URL, this sends them index.html
@app.route('/')
def index():
    return render_template('index.html')


# ── Route 2: Receive uploaded files & run analysis ────────────────
# When the user clicks Analyze, the browser sends the files here
@app.route('/analyze', methods=['POST'])
def analyze():
    try:
        return _analyze_inner()
    except Exception as e:
        import traceback
        return jsonify({
            'error': 'Unexpected server error',
            'details': traceback.format_exc()
        }), 500


def _analyze_inner():
    # ── Step 1: Get the uploaded files ──
    if 'alphafold' not in request.files or 'bound' not in request.files:
        return jsonify({'error': 'Please upload both PDB files.'}), 400

    alphafold_file = request.files['alphafold']
    bound_file     = request.files['bound']
    ligand_code    = request.form.get('ligand_code', 'IRE').strip().upper()
    protein_name = request.form.get('protein_name', 'Unknown Protein').strip()
    mutation_name = request.form.get('mutation_name', 'Unknown Mutation').strip()
    disease_name = request.form.get('disease_name', 'Unknown Disease').strip()

    # ── Step 2: Save with unique names so simultaneous users don't overwrite each other ──
    session_id     = str(uuid.uuid4())[:8]   # random 8-character ID
    alphafold_path = os.path.join(UPLOAD_FOLDER, f"{session_id}_alphafold.pdb")
    bound_path     = os.path.join(UPLOAD_FOLDER, f"{session_id}_bound.pdb")
    results_dir    = os.path.join(RESULTS_FOLDER, session_id)
    images_dir     = os.path.join(IMAGES_FOLDER,  session_id)

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(images_dir,  exist_ok=True)

    alphafold_file.save(alphafold_path)
    bound_file.save(bound_path)
    # Save analysis to history
    history_path = os.path.join(BASE_DIR, "data", "history.csv")
    history_exists = os.path.exists(history_path)

    with open(history_path, "a", newline="", encoding='utf-8', errors='replace') as f:
        writer = csv.writer(f)

        if not history_exists:
            writer.writerow([
                "session_id",
                "protein_name",
                "mutation_name",
                "disease_name",
                "ligand_code",
                "alphafold_file",
                "bound_file",
                "timestamp"
            ])

        writer.writerow([
            session_id,
            protein_name,
            mutation_name,
            disease_name,
            ligand_code,
            alphafold_file.filename,
            bound_file.filename,
            datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ])

    python = sys.executable  # use the same Python that's running Flask

    # ── Step 3: Run analyze_structure.py ──
    try:
        subprocess.run([
            python,
            os.path.join(SCRIPTS_FOLDER, 'analyze_structure.py'),
            alphafold_path,
            bound_path,
            results_dir,
            ligand_code
        ], check=True, capture_output=True, text=True)

        subprocess.run([
            python,
            os.path.join(SCRIPTS_FOLDER, 'find_binding_pocket.py'),
            bound_path,
            ligand_code,
            results_dir
        ], check=True, capture_output=True, text=True)

        subprocess.run([
            python,
            os.path.join(SCRIPTS_FOLDER, 'plot_summary.py'),
            results_dir,
            images_dir
        ], check=True, capture_output=True, text=True)

    except subprocess.CalledProcessError as e:
        return jsonify({
            'error': 'Backend script failed',
            'details': e.stderr
        }), 500

    # ── Step 4: Run new analysis scripts (non-fatal — each wrapped individually) ──
    new_scripts = [
        [python, os.path.join(SCRIPTS_FOLDER, 'analyze_hbonds.py'),
         bound_path, ligand_code, results_dir],
        [python, os.path.join(SCRIPTS_FOLDER, 'analyze_occupancy.py'),
         bound_path, ligand_code, results_dir],
        [python, os.path.join(SCRIPTS_FOLDER, 'analyze_metal_coordination.py'),
         bound_path, ligand_code, results_dir],
        [python, os.path.join(SCRIPTS_FOLDER, 'ramachandran_plot.py'),
         alphafold_path, bound_path, images_dir, results_dir],
        [python, os.path.join(SCRIPTS_FOLDER, 'sequence_comparison.py'),
         alphafold_path, bound_path, results_dir, ligand_code],
    ]
    for cmd in new_scripts:
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=120)
        except Exception:
            pass  # non-fatal; partial results still returned

    # ── Step 5: Read the results and send them back to the browser ──

    # Read the CSV stats
    csv_path = os.path.join(results_dir, "structure_summary.csv")
    stats = []
    if os.path.exists(csv_path):
        with open(csv_path, encoding='utf-8', errors='replace') as f:
            reader = csv.DictReader(f)
            stats = list(reader)
            
    advanced_path = os.path.join(results_dir, "advanced_metrics.csv")
    advanced_metrics = []
    if os.path.exists(advanced_path):
        with open(advanced_path, encoding='utf-8', errors='replace') as f:
            reader = csv.DictReader(f)
            advanced_metrics = list(reader)

    # Read the binding pocket residues text file
    pocket_path = os.path.join(results_dir, "binding_pocket_residues.txt")
    pocket_text = ""
    if os.path.exists(pocket_path):
        with open(pocket_path, encoding='utf-8', errors='replace') as f:
            pocket_text = f.read()

    # Convert the graph PNG to base64 (a way to embed an image directly in JSON)
    graphs = []
    for filename in os.listdir(images_dir):
        if filename.endswith('.png'):
            filepath = os.path.join(images_dir, filename)
            with open(filepath, 'rb') as f:
                encoded = base64.b64encode(f.read()).decode('utf-8')
                graphs.append({'name': filename, 'data': encoded})
    # Read ML summary report
    ml_summary_path = os.path.join(RESULTS_FOLDER, "ml_summary_report.txt")
    ml_summary_text = ""

    if os.path.exists(ml_summary_path):
        with open(ml_summary_path, "r", encoding='utf-8', errors='replace') as f:
            ml_summary_text = f.read()


    # Read global ML images
    ml_graphs = []

    global_ml_images = [
        "pca_clusters_improved.png",
        "mutation_similarity_heatmap.png",
        "feature_importance.png",
        "algorithm_comparison.png",
        "feature_correlation_matrix.png"
    ]

    for filename in global_ml_images:
        filepath = os.path.join(IMAGES_FOLDER, filename)

        if os.path.exists(filepath):
            with open(filepath, "rb") as f:
                encoded = base64.b64encode(f.read()).decode("utf-8")
                ml_graphs.append({
                    "name": filename,
                    "data": encoded
                })
    # ── Read new analysis results ──

    def read_csv_safe(path):
        if not os.path.exists(path):
            return []
        with open(path, newline='', encoding='utf-8', errors='replace') as f:
            return list(csv.DictReader(f))

    def read_txt_safe(path):
        if not os.path.exists(path):
            return ''
        with open(path, encoding='utf-8', errors='replace') as f:
            return f.read()

    hbond_data       = read_csv_safe(os.path.join(results_dir, 'hbond_interactions.csv'))
    hbond_summary    = read_txt_safe(os.path.join(results_dir, 'hbond_summary.txt'))
    occ_low          = read_csv_safe(os.path.join(results_dir, 'low_occupancy_residues.csv'))
    occ_high         = read_csv_safe(os.path.join(results_dir, 'high_bfactor_residues.csv'))
    occ_summary      = read_txt_safe(os.path.join(results_dir, 'occupancy_bfactor_summary.txt'))
    metal_data       = read_csv_safe(os.path.join(results_dir, 'metal_coordination.csv'))
    metal_summary    = read_txt_safe(os.path.join(results_dir, 'metal_coordination_summary.txt'))
    seq_comp_data    = read_csv_safe(os.path.join(results_dir, 'sequence_comparison.csv'))
    seq_comp_summary = read_txt_safe(os.path.join(results_dir, 'sequence_comparison_summary.txt'))
    res_corr_data    = read_csv_safe(os.path.join(results_dir, 'residue_correspondence.csv'))

    # Send everything back to the browser as JSON
    return jsonify({
        'status': 'success',
        'stats': stats,
        'advanced_metrics': advanced_metrics,
        'pocket': pocket_text,
        'graphs': graphs,
        'ml_summary': ml_summary_text,
        'ml_graphs': ml_graphs,
        'hbond_data': hbond_data,
        'hbond_summary': hbond_summary,
        'occupancy_low': occ_low,
        'occupancy_high': occ_high,
        'occupancy_summary': occ_summary,
        'metal_data': metal_data,
        'metal_summary': metal_summary,
        'seq_comparison': seq_comp_data,
        'seq_comparison_summary': seq_comp_summary,
        'residue_correspondence': res_corr_data,
        'alphafold_id': f"{session_id}_alphafold",
        'bound_id': f"{session_id}_bound",
        'results_folder': session_id,
        'metadata': {
            'protein_name': protein_name,
            'mutation_name': mutation_name,
            'ligand_code': ligand_code,
            'disease_name': disease_name
        },
    })


# ── Route 3: Serve uploaded PDB files to the 3D viewer ────────────
# The NGL viewer in the browser needs to fetch the PDB file as a URL
@app.route('/pdb/<filename>')
def serve_pdb(filename):
    from flask import send_from_directory
    return send_from_directory(UPLOAD_FOLDER, filename + '.pdb',
                               mimetype='text/plain')

@app.route('/download/<session_id>/<filename>')
def download_file(session_id, filename):

    from flask import send_from_directory

    folder = os.path.join(RESULTS_FOLDER, session_id)

    return send_from_directory(folder, filename, as_attachment=True)
# ── Start the server ──────────────────────────────────────────────
@app.route('/history')
def history():
    history_path = os.path.join(BASE_DIR, "data", "history.csv")

    rows = []

    if os.path.exists(history_path):
        with open(history_path, encoding='utf-8', errors='replace') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

    rows.reverse()

    return jsonify(rows)

@app.route('/blast', methods=['POST'])
def blast_search():
    try:
        return _blast_inner()
    except Exception as e:
        import traceback
        return jsonify({'error': 'BLAST server error', 'details': traceback.format_exc()}), 500


def _blast_inner():
    """Separate endpoint for NCBI BLAST — called by the frontend after main analysis."""
    session_id  = request.form.get('session_id', '').strip()
    pdb_file    = request.form.get('pdb_file',   '').strip()
    threshold   = request.form.get('threshold',  '80').strip()

    if not session_id or not pdb_file:
        return jsonify({'error': 'Missing session_id or pdb_file'}), 400

    bound_path  = os.path.join(UPLOAD_FOLDER,  pdb_file + '.pdb')
    results_dir = os.path.join(RESULTS_FOLDER, session_id)

    if not os.path.exists(bound_path):
        return jsonify({'error': 'PDB file not found on server'}), 404

    python = sys.executable
    try:
        subprocess.run([
            python,
            os.path.join(SCRIPTS_FOLDER, 'blast_search.py'),
            bound_path, results_dir, threshold,
        ], check=True, capture_output=True, text=True, timeout=180)
    except subprocess.TimeoutExpired:
        return jsonify({'error': 'BLAST timed out (>180 s). Try a shorter sequence.'}), 504
    except subprocess.CalledProcessError as e:
        return jsonify({'error': 'BLAST failed', 'details': e.stderr}), 500

    blast_csv = os.path.join(results_dir, 'blast_results.csv')
    blast_txt = os.path.join(results_dir, 'blast_summary.txt')

    blast_data = []
    if os.path.exists(blast_csv):
        with open(blast_csv, newline='', encoding='utf-8', errors='replace') as f:
            blast_data = list(csv.DictReader(f))

    blast_summary = ''
    if os.path.exists(blast_txt):
        with open(blast_txt, encoding='utf-8', errors='replace') as f:
            blast_summary = f.read()

    return jsonify({
        'status': 'success',
        'blast_data': blast_data,
        'blast_summary': blast_summary,
    })


if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get('PORT', 5000)))
